## P1 readings

- designs 23400: matched 4317, unmatched 19083; worst probe−engine design difference 5.55E-016
- (i) unmatched with rho > 1e-13: 0; of them cleared by the zeroed recheck (rho ≤ 1e-13): 0; max rho among unmatched 2.24E-017
- (ii) matched rho: min 1.53E-017, 1% 5.19E-008, median 3.38E-002; below 1e-12: 2, below 1e-10: 6

| matching | engine verdict | designs | of them N > n |
|---|---|---|---|
| false | ThinSupport | 19083 | 10531 |
| true | None | 787 | 787 |
| true | ThinSupport | 3530 | 2722 |

Matched, N > n, refused by the engine, by reason and rho band:

- ThinSupport: 2722; rho < 1e-10: 4, 1e-10..1e-6: 49, ≥ 1e-6: 2669; by site law uniform 1167, clustered 418, gapped 1137

### Products (6000 designs)

- full rank with Mx or My false (per-axis necessity): 0
- Mx and My true but rank-deficient (per-axis sufficiency fails): 3040 of 4211 with both matched
- full rank without a bipartite matching (its necessity): 0
- bipartite matching but rank-deficient (its sufficiency fails): 378 of 1549 with a matching

| sites | designs | full rank | Mx∧My | bipartite | Mx∧My ∧ deficient | bipartite ∧ deficient |
|---|---|---|---|---|---|---|
| uniform | 1500 | 681 | 1335 | 681 | 654 | 0 |
| diagonal | 1500 | 0 | 1333 | 374 | 1333 | 374 |
| near-diagonal | 1500 | 363 | 1353 | 366 | 990 | 3 |
| clustered | 1500 | 127 | 190 | 128 | 63 | 1 |

Kept counterexamples (bipartite matching, rank-deficient; first 8):

- 0,0,16,16,diagonal,0,1,1,1,1.576E-020
- 0,0,16,16,diagonal,1,1,1,1,1.693E-018
- 0,0,16,16,diagonal,2,1,1,1,1.708E-018
- 0,0,16,16,diagonal,3,1,1,1,8.666E-019
- 0,0,16,16,diagonal,4,1,1,1,8.317E-019
- 0,0,16,16,diagonal,5,1,1,1,6.887E-019
- 0,0,16,16,diagonal,6,1,1,1,2.082E-019
- 0,0,16,16,diagonal,7,1,1,1,3.821E-019
