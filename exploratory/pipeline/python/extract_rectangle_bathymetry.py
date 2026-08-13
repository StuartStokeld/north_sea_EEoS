#!/usr/bin/env python3
"""Zonal GEBCO bathymetry + TID QC for H2/H3 rectangles.

Writes a CSV with mean depth, depth range, local depth-gradient direction/
magnitude at the centroid, percent sounding coverage, and distance to the
GEBCO grid edge.

Run from repo root (or via the R orchestrator):
  python3 pipeline/python/extract_rectangle_bathymetry.py \\
    --shapefile gis/ICES_rectangles/ICES_Statistical_Rectangles_Eco.shp \\
    --panel-ids outputs/_tmp_panel_stat_rec.txt \\
    --bathy data/external/.../gebco_2026_..._geotiff.tif \\
    --tid data/external/.../gebco_2026_tid_..._geotiff.tif \\
    --out outputs/bathymetry_by_rectangle.csv
"""

from __future__ import annotations

import argparse
import csv
import math
import sys

import numpy as np
from osgeo import gdal, ogr, osr

gdal.UseExceptions()
ogr.UseExceptions()

# GEBCO TID: direct measurements (Table 2)
DIRECT_TID = {10, 11, 12, 13, 14, 15, 16, 17}
LAND_TID = {0}

# Window half-width (cells) for centroid gradient (~21 x 21 at 15")
GRAD_HALF = 10


def read_ids(path: str) -> list[str]:
    with open(path, encoding="utf-8") as f:
        return [ln.strip() for ln in f if ln.strip()]


def open_rectangles(shapefile: str, ids: set[str]):
    ds = ogr.Open(shapefile)
    if ds is None:
        raise RuntimeError(f"Cannot open shapefile: {shapefile}")
    layer = ds.GetLayer(0)
    feats = []
    for feat in layer:
        name = str(feat.GetField("ICESNAME")).strip()
        if name in ids:
            geom = feat.GetGeometryRef().Clone()
            feats.append((name, geom))
    if len(feats) != len(ids):
        missing = sorted(ids - {n for n, _ in feats})
        raise RuntimeError(
            f"Shapefile matched {len(feats)}/{len(ids)} ids; missing e.g. {missing[:10]}"
        )
    # Preserve caller order later via dict
    return {n: g for n, g in feats}, ds  # keep ds alive


def raster_window(ds, minx, maxx, miny, maxy, pad_cells: int = 0):
    gt = ds.GetGeoTransform()
    inv = gdal.InvGeoTransform(gt)
    # pixel coords of corners
    cols, rows = [], []
    for x, y in ((minx, maxy), (maxx, maxy), (minx, miny), (maxx, miny)):
        px, py = gdal.ApplyGeoTransform(inv, x, y)
        cols.append(px)
        rows.append(py)
    xoff = int(math.floor(min(cols))) - pad_cells
    yoff = int(math.floor(min(rows))) - pad_cells
    x2 = int(math.ceil(max(cols))) + pad_cells
    y2 = int(math.ceil(max(rows))) + pad_cells
    xoff = max(0, xoff)
    yoff = max(0, yoff)
    x2 = min(ds.RasterXSize, x2)
    y2 = min(ds.RasterYSize, y2)
    win_xsize = max(0, x2 - xoff)
    win_ysize = max(0, y2 - yoff)
    return xoff, yoff, win_xsize, win_ysize


def pixel_mask_for_geom(geom, gt, xoff, yoff, win_xsize, win_ysize):
    """Rasterize polygon into a boolean mask aligned with the window."""
    if win_xsize == 0 or win_ysize == 0:
        return np.zeros((0, 0), dtype=bool)
    mem = gdal.GetDriverByName("MEM").Create("", win_xsize, win_ysize, 1, gdal.GDT_Byte)
    # Window geotransform
    win_gt = (
        gt[0] + xoff * gt[1] + yoff * gt[2],
        gt[1],
        gt[2],
        gt[3] + xoff * gt[4] + yoff * gt[5],
        gt[4],
        gt[5],
    )
    mem.SetGeoTransform(win_gt)
    # burn
    drv = ogr.GetDriverByName("MEM")
    tmp = drv.CreateDataSource("")
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    lyr = tmp.CreateLayer("poly", srs=srs, geom_type=ogr.wkbPolygon)
    fdefn = lyr.GetLayerDefn()
    feat = ogr.Feature(fdefn)
    feat.SetGeometry(geom)
    lyr.CreateFeature(feat)
    gdal.RasterizeLayer(mem, [1], lyr, burn_values=[1])
    mask = mem.GetRasterBand(1).ReadAsArray().astype(bool)
    return mask


def window_geotransform(gt, xoff, yoff):
    return (
        gt[0] + xoff * gt[1] + yoff * gt[2],
        gt[1],
        gt[2],
        gt[3] + xoff * gt[4] + yoff * gt[5],
        gt[4],
        gt[5],
    )


def sample_gradient_at_centroid(bathy_ds, lon, lat):
    """Local depth gradient at centroid.

    Depth = -elevation (positive deeper). Gradient uses metres-per-degree
    with longitude compressed by cos(lat).

    Returns (bearing_cross_deg, bearing_along_deg, magnitude, elev_at_centroid)
    Bearings use gstat plane convention: degrees counterclockwise from +x (east),
    so 0=east, 90=north.
    """
    gt = bathy_ds.GetGeoTransform()
    inv = gdal.InvGeoTransform(gt)
    px, py = gdal.ApplyGeoTransform(inv, lon, lat)
    cx, cy = int(round(px)), int(round(py))
    xoff = max(0, cx - GRAD_HALF)
    yoff = max(0, cy - GRAD_HALF)
    x2 = min(bathy_ds.RasterXSize, cx + GRAD_HALF + 1)
    y2 = min(bathy_ds.RasterYSize, cy + GRAD_HALF + 1)
    arr = bathy_ds.ReadAsArray(xoff, yoff, x2 - xoff, y2 - yoff).astype(float)
    if arr.size < 9:
        return (np.nan, np.nan, np.nan, np.nan)

    depth = -arr  # positive deeper
    # cell size in degrees
    dlon = abs(gt[1])
    dlat = abs(gt[5])
    # metres approx
    lat_rad = math.radians(lat)
    m_per_deg_lat = 111320.0
    m_per_deg_lon = 111320.0 * math.cos(lat_rad)
    # np.gradient: axis0 = rows (north->south if gt[5]<0), axis1 = cols (west->east)
    # With north-up GeoTIFF, row increase = south = decreasing lat
    g_row, g_col = np.gradient(depth, dlat * m_per_deg_lat, dlon * m_per_deg_lon)
    # Convert row gradient (per metre south) to per-metre north: negate
    ddepth_dy = -g_row  # north component
    ddepth_dx = g_col  # east component

    ic = cy - yoff
    jc = cx - xoff
    ic = min(max(ic, 0), ddepth_dy.shape[0] - 1)
    jc = min(max(jc, 0), ddepth_dy.shape[1] - 1)
    # smooth local mean over 3x3 around centre
    r0, r1 = max(0, ic - 1), min(ddepth_dy.shape[0], ic + 2)
    c0, c1 = max(0, jc - 1), min(ddepth_dy.shape[1], jc + 2)
    dx = float(np.nanmean(ddepth_dx[r0:r1, c0:c1]))
    dy = float(np.nanmean(ddepth_dy[r0:r1, c0:c1]))
    mag = math.hypot(dx, dy)
    # angle from +east toward +north (counterclockwise): atan2(north, east)
    bearing_cross = math.degrees(math.atan2(dy, dx)) % 360.0
    bearing_along = (bearing_cross + 90.0) % 360.0
    elev = float(arr[ic, jc])
    return bearing_cross, bearing_along, mag, elev


def dist_to_gebco_edge(lon, lat, xmin, xmax, ymin, ymax):
    return float(min(lon - xmin, xmax - lon, lat - ymin, ymax - lat))


def process(args):
    ids_list = read_ids(args.panel_ids)
    ids = set(ids_list)
    feats, _shp_ds = open_rectangles(args.shapefile, ids)

    bathy = gdal.Open(args.bathy)
    tid = gdal.Open(args.tid)
    gt = bathy.GetGeoTransform()
    # GEBCO extent from geotransform + size
    xmin = gt[0]
    ymax = gt[3]
    xmax = gt[0] + bathy.RasterXSize * gt[1]
    ymin = gt[3] + bathy.RasterYSize * gt[5]
    # ensure ordering
    xmin, xmax = min(xmin, xmax), max(xmin, xmax)
    ymin, ymax = min(ymin, ymax), max(ymin, ymax)

    rows_out = []
    for name in ids_list:
        geom = feats[name]
        env = geom.GetEnvelope()  # minx, maxx, miny, maxy
        minx, maxx, miny, maxy = env
        lon = 0.5 * (minx + maxx)
        lat = 0.5 * (miny + maxy)

        xoff, yoff, wx, wy = raster_window(bathy, minx, maxx, miny, maxy, pad_cells=1)
        elev = bathy.ReadAsArray(xoff, yoff, wx, wy).astype(float)
        tid_arr = tid.ReadAsArray(xoff, yoff, wx, wy)
        mask = pixel_mask_for_geom(geom, gt, xoff, yoff, wx, wy)
        if mask.shape != elev.shape:
            # reshape safety
            mask = mask[: elev.shape[0], : elev.shape[1]]

        cells = elev[mask]
        tid_cells = tid_arr[mask]
        if cells.size == 0:
            mean_elev = np.nan
            depth_min = np.nan
            depth_max = np.nan
            depth_range = np.nan
            mean_depth = np.nan
            pct_sounding = np.nan
            n_cells = 0
            n_nonland = 0
            n_direct = 0
        else:
            mean_elev = float(np.mean(cells))
            # depth positive below sea level
            depth = -cells
            depth_min = float(np.min(depth))
            depth_max = float(np.max(depth))
            depth_range = depth_max - depth_min
            mean_depth = float(np.mean(depth))
            n_cells = int(cells.size)
            nonland = ~np.isin(tid_cells, list(LAND_TID))
            n_nonland = int(np.sum(nonland))
            n_direct = int(np.sum(np.isin(tid_cells, list(DIRECT_TID))))
            pct_sounding = (
                100.0 * n_direct / n_nonland if n_nonland > 0 else np.nan
            )

        b_cross, b_along, gmag, elev_c = sample_gradient_at_centroid(bathy, lon, lat)
        edge = dist_to_gebco_edge(lon, lat, xmin, xmax, ymin, ymax)

        rows_out.append(
            {
                "stat_rec": name,
                "lon": lon,
                "lat": lat,
                "mean_elevation_m": mean_elev,
                "mean_depth_m": mean_depth,
                "depth_min_m": depth_min,
                "depth_max_m": depth_max,
                "depth_range_m": depth_range,
                "grad_bearing_cross_deg": b_cross,
                "grad_bearing_along_deg": b_along,
                "grad_magnitude_m_per_m": gmag,
                "elev_at_centroid_m": elev_c,
                "pct_sounding": pct_sounding,
                "n_cells": n_cells,
                "n_nonland": n_nonland,
                "n_direct": n_direct,
                "dist_to_gebco_edge_deg": edge,
            }
        )

    fieldnames = list(rows_out[0].keys())
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows_out:
            w.writerow(r)
    print(f"Wrote {len(rows_out)} rows to {args.out}", file=sys.stderr)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--shapefile", required=True)
    p.add_argument("--panel-ids", required=True)
    p.add_argument("--bathy", required=True)
    p.add_argument("--tid", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()
    process(args)


if __name__ == "__main__":
    main()
