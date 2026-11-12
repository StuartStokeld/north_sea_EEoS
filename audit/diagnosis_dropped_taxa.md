# Investigation: taxa dropped for lack of length–weight coefficients (F-B8)

**Status:** investigation only; nothing changed in the pipeline. Scripts:

- `audit/scripts/rebuild_lw_lookup.R` (draft, re-runnable)
- `audit/scripts/dropped_taxa_effect.R`

Outputs are in `audit/outputs/lw_rebuild/` and `audit/outputs/dropped_taxa/`.

## The dropped taxa

`add_lw_mass_to_bins()` inner-joins the HL bins to `outputs/fishbase_lw_lookup_v2.csv` on the DATRAS AphiaID.
Any taxon without finite a and b drops out of S, N, E and m_min.

In the 10,128 H1 hauls, **53 taxa** are dropped:

- 29 are absent from the lookup altogether. I counted 30 at GATE A over all effort-kept hauls; 29 occur in H1
  hauls.
- 24 are in the lookup with NA a/b. These are mostly the "32 MISSING names" in the methods.

Together they account for **0.39% of individuals**. **4,902 hauls lose at least one species**: median 8.3% of
the haul's species, maximum 33%.

**Why they're absent.** The lookup's species list (`all_spp` in the notebook) came from FishGlob's accepted
names. The AphiaIDs that FishGlob lumped or renamed never got an entry under the DATRAS AphiaID the pipeline
joins on.

Largest by individuals (full list with WoRMS status in `dropped_taxa/dropped_taxa.csv`):

| AphiaID | Taxon (WoRMS) | Rank | Individuals | Hauls | Resolved by rebuild |
|---|---|---|---|---|---|
| 126716 | *Argentina sphyraena* | species | 77,518 | 1,843 | TL record |
| 127387 | *Syngnathus acus* | species | 75,951 | 256 | TL record |
| 126793 | *Callionymus maculatus* | species | 25,146 | 1,652 | TL record |
| 125537 | Gobiidae | family | 21,686 | 460 | family mean |
| 125516 | Ammodytidae | family | 15,545 | 258 | family mean |
| 126715 | *Argentina silus* | species | 7,992 | 343 | TL record |
| 126928 | *Pomatoschistus minutus* | species | 6,025 | 291 | TL record |
| 126415 | *Alosa fallax* | species | 3,546 | 344 | TL record |
| 127389 | *Syngnathus rostellatus* | species | 3,418 | 264 | TL record |
| 127219 | *Liparis liparis* | species | 2,463 | 372 | TL record |
| 127147 | *Phrynorhombus norvegicus* | species | 942 | 136 | TL record |

## Name matching

All 53 AphiaIDs are **accepted** in WoRMS, except *Dipturus flossada* (taxon inquirendum). WoRMS was queried for
all 212 fish AphiaIDs in the pipeline's HL rows; the records are cached in `lw_rebuild/worms_records.csv`.
WoRMS returned HTTP 503 for part of the session, so the script falls back to DATRAS species names and caches
only successful records.

WoRMS valid names matched FishBase valid names directly for every species-rank taxon that resolved. The
**FishBase synonym step was dropped by scope decision** and is logged as a follow-up. It matters for names like
*Phrynorhombus norvegicus*, which FishGlob's table holds as *Zeugopterus norvegicus*.

## Rebuilt lookup (FishBase release 26.06, rfishbase 5.0.3)

The rebuild uses the notebook's rules, replicated: TL > FL > SL > untyped > other, then the largest `Number`;
then genus mean, then family mean over TL and untyped records of the listed species. It adds a FishBase-wide
genus/family tier for anything still unresolved.

Results for the 212 taxa:

| Resolution | Taxa |
|---|---|
| TL record | 157 |
| Genus mean | 25 |
| Family mean | 20 |
| Genus / family mean (FishBase-wide tier) | 1 / 1 |
| SL record (needs conversion) | 4 |
| Untyped / other | 1 / 1 |
| Missing | 2 (*Ostrea edulis*, *Todarodes*: not fish, see F-B12) |

So **51 of the 53 dropped taxa resolve.**

**Comparison with the current lookup.** 157 of the 177 taxa with coefficients are present in the rebuild; the
other 20 are keyed on FishGlob AphiaIDs that never occur in the pipeline's HL. Of those 157:

- 110 have identical a and b.
- 25 change the mass of a 30 cm fish by more than 10%. The causes are FishBase release differences (current
  lookup: "accessed May 2026") and ties in record choice.

Release 25.04 couldn't be downloaded (S3 timeout), so the exact original release is unconfirmed.

## Effect on S, N, E and H1 (EEoS re-solved for changed hauls)

| Variant | Hauls changed | S | Median N ratio | Median E ratio | Failed solves | 1:1 log R² | Bias | Bias-removed R² | Median B_pred/B_obs | cor² |
|---|---|---|---|---|---|---|---|---|---|---|
| Current lookup | — | — | — | — | 5 | 0.6610 | −0.652 | 0.9282 | 1.958 | 0.9713 |
| Current + 51 resolved dropped taxa | 4,902 | +1 (median) in 4,902 | 1.0023 | 1.0026 | +2 | 0.6580 | −0.657 | 0.9287 | 1.965 | 0.9710 |
| Fully rebuilt lookup | 9,343 | +1 in 4,902 | 1.0001 | 1.0002 | +3 | 0.6584 | −0.656 | 0.9283 | 1.963 | 0.9710 |

**Reading.** The H1 effect is small. The main consequence is that S is understated by one species in about half
of all hauls, which matters for the METE inputs and for any analysis of S.

H2/H3 were not refitted here. That is deferred to the fix branch, where the lookup change would be combined with
F-A1/A2.
