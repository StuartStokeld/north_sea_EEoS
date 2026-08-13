# Spec A+B — VIF diagnostic

Rectangle-year grain (unique `stat_rec` × `year` in analysis data).
VIF = 1/(1−R²) from each term ~ remaining continuous terms.

| Term | VIF | R² (others) | n |
|------|-----|-------------|---|
| FP_between | 3.131 | 0.6806 | 4768 |
| FP_within | 1.026 | 0.0249 | 4768 |
| FP_between_lag | 3.641 | 0.7254 | 4768 |
| B_lag_neighbour | 1.457 | 0.3136 | 4768 |

