### Graphs

| Graph | Arm | Filter | Edges (K̂) | Largest component before repair | Repair edges | Components | Bandwidth rule | a | Cut by the floor | Scale edges | Gauge |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Iris, K = 5, Gaussian | Base | mutual | 262 (3.49) | 53% | 15 | 1 | mean edge | 0.3065 | 0 | after | 0.2863 |
|  | OR filter | OR | 517 (6.89) | 65% | 1 | 1 | mean edge | 0.3896 | 0 | after | 0.1451 |
|  | No repair | mutual | 247 (3.29) | 53% | 0 | 16 | mean edge | 0.2893 | 0 | after | 0.3036 |
|  | MAD bandwidth | mutual | 262 (3.49) | 53% | 15 | 1 | MAD | 0.1212 | 0 | after | 0.2863 |
|  | Scale before repair | mutual | 262 (3.49) | 53% | 15 | 1 | mean edge | 0.2893 | 0 | before | 0.2863 |
|  | Compiler defaults | OR | 504 (6.72) | 65% | 1 | 1 | MAD | 0.1212 | 13 | before | none |
| Iris, K = 10, Gaussian | Base | mutual | 528 (7.04) | 63% | 5 | 1 | mean edge | 0.3645 | 0 | after | 0.1421 |
|  | OR filter | OR | 995 (13.27) | 67% | 1 | 1 | mean edge | 0.4754 | 0 | after | 0.0754 |
|  | No repair | mutual | 523 (6.97) | 63% | 0 | 6 | mean edge | 0.3602 | 0 | after | 0.1434 |
|  | MAD bandwidth | mutual | 525 (7.00) | 63% | 5 | 1 | MAD | 0.1414 | 3 | after | 0.1429 |
|  | Scale before repair | mutual | 528 (7.04) | 63% | 5 | 1 | mean edge | 0.3602 | 0 | before | 0.1421 |
|  | Compiler defaults | OR | 967 (12.89) | 67% | 1 | 1 | MAD | 0.1414 | 28 | before | none |
| Iris, K = 10, Cauchy | Base | mutual | 528 (7.04) | 63% | 5 | 1 | mean edge | 0.3645 | 0 | after | 0.1421 |
|  | OR filter | OR | 995 (13.27) | 67% | 1 | 1 | mean edge | 0.4754 | 0 | after | 0.0754 |
|  | No repair | mutual | 523 (6.97) | 63% | 0 | 6 | mean edge | 0.3602 | 0 | after | 0.1434 |
|  | MAD bandwidth | mutual | 528 (7.04) | 63% | 5 | 1 | MAD | 0.09535 | 0 | after | 0.1421 |
|  | Scale before repair | mutual | 528 (7.04) | 63% | 5 | 1 | mean edge | 0.3602 | 0 | before | 0.1421 |
|  | Compiler defaults | OR | 995 (13.27) | 67% | 1 | 1 | MAD | 0.09535 | 0 | before | none |
| The toy, K = 10, Gaussian | Base | mutual | 12,655 (7.91) | 99% | 12 | 1 | mean edge | 0.1622 | 0 | after | 0.1264 |
|  | OR filter | OR | 19,357 (12.10) | 100% | 0 | 1 | mean edge | 0.2206 | 0 | after | 0.0827 |
|  | No repair | mutual | 12,643 (7.90) | 99% | 0 | 13 | mean edge | 0.162 | 0 | after | 0.1265 |
|  | MAD bandwidth | mutual | 12,019 (7.51) | 99% | 12 | 56 | MAD | 0.05568 | 636 | after | 0.1331 |
|  | Scale before repair | mutual | 12,655 (7.91) | 99% | 12 | 1 | mean edge | 0.162 | 0 | before | 0.1264 |
|  | Compiler defaults | OR | 14,747 (9.22) | 100% | 0 | 57 | MAD | 0.05568 | 2109 | before | none |
| Landsat, K = 5, Gaussian | Base | mutual | 9,048 (2.81) | 78% | 983 | 1 | mean edge | 22.71 | 0 | after | 0.3556 |
|  | OR filter | OR | 24,284 (7.55) | 99% | 0 | 1 | mean edge | 24.8 | 0 | after | 0.1325 |
|  | No repair | mutual | 8,065 (2.51) | 78% | 0 | 984 | mean edge | 22.32 | 0 | after | 0.3990 |
|  | MAD bandwidth | mutual | 8,925 (2.77) | 78% | 983 | 66 | MAD | 6.732 | 123 | after | 0.3605 |
|  | Scale before repair | mutual | 9,048 (2.81) | 78% | 983 | 1 | mean edge | 22.32 | 0 | before | 0.3556 |
|  | Compiler defaults | OR | 15,742 (4.89) | 99% | 0 | 177 | MAD | 6.732 | 524 | before | none |
| Landsat, K = 10, Gaussian | Base | mutual | 17,109 (5.32) | 95% | 282 | 1 | mean edge | 23.84 | 0 | after | 0.1881 |
|  | OR filter | OR | 47,850 (14.87) | 100% | 0 | 1 | mean edge | 26.59 | 0 | after | 0.0672 |
|  | No repair | mutual | 16,827 (5.23) | 95% | 0 | 283 | mean edge | 23.77 | 0 | after | 0.1912 |
|  | MAD bandwidth | mutual | 16,872 (5.24) | 95% | 282 | 36 | MAD | 7.762 | 237 | after | 0.1907 |
|  | Scale before repair | mutual | 17,109 (5.32) | 95% | 282 | 1 | mean edge | 23.77 | 0 | before | 0.1881 |
|  | Compiler defaults | OR | 31,459 (9.78) | 100% | 0 | 51 | MAD | 7.762 | 853 | before | none |
| BlattHierarchy, K = 5, Cauchy | Base | mutual | 2,169 (3.62) | 4% | 82 | 1 | mean edge | 0.08313 | 0 | after | 0.2766 |
|  | OR filter | OR | 3,933 (6.55) | 8% | 20 | 1 | mean edge | 0.1059 | 0 | after | 0.1526 |
|  | No repair | mutual | 2,087 (3.48) | 4% | 0 | 83 | mean edge | 0.06391 | 0 | after | 0.2875 |
|  | MAD bandwidth | mutual | 2,169 (3.62) | 4% | 82 | 1 | MAD | 0.02741 | 0 | after | 0.2766 |
|  | Scale before repair | mutual | 2,169 (3.62) | 4% | 82 | 1 | mean edge | 0.06391 | 0 | before | 0.2766 |
|  | Compiler defaults | OR | 3,933 (6.55) | 8% | 20 | 1 | MAD | 0.02741 | 0 | before | none |

### auto

| Graph | Base | OR filter | No repair | Scale before repair | MAD bandwidth | No gauge | q = 5 | q = 40 | Compiler defaults |
|---|---|---|---|---|---|---|---|---|---|
| Iris, K = 5, Gaussian | T_ps [0.00446, 0.357] | T_ps [0.00446, 0.357] | T_ps [0.00446, 0.357] | T_ps [0.00446, 0.357] | [1.16·10⁻⁴, 0.0116] | [0.0222, 2.22] | T_ps [0.00646, 0.516] | T_ps [0.00381, 0.305] | [1.2·10⁻⁴, 0.012] |
| Iris, K = 10, Gaussian | T_ps [0.00446, 0.357] | T_ps [0.00446, 0.357] | T_ps [0.00446, 0.357] | T_ps [0.00446, 0.357] | [1.1·10⁻⁴, 0.011] | [0.0442, 4.42] | T_ps [0.00646, 0.516] | T_ps [0.00381, 0.305] | [1.41·10⁻⁴, 0.0141] |
| Iris, K = 10, Cauchy | [0.00451, 0.451] | [0.00477, 0.477] | [0.00461, 0.461] | [0.00441, 0.441] | [1.11·10⁻⁴, 0.0111] | [0.0317, 3.17] | [0.00839, 0.839] | [0.00366, 0.366] | [6.62·10⁻⁴, 0.0662] |
| The toy, K = 10, Gaussian | T_ps [0.00446, 0.357] | T_ps [0.00446, 0.357] | T_ps [0.00446, 0.357] | T_ps [0.00446, 0.357] | [1.05·10⁻⁴, 0.0105] | [0.0564, 5.64] | T_ps [0.00646, 0.516] | T_ps [0.00381, 0.305] | [3.58·10⁻⁴, 0.0358] |
| Landsat, K = 5, Gaussian | T_ps [0.00446, 0.357] | T_ps [0.00446, 0.357] | T_ps [0.00446, 0.357] | T_ps [0.00446, 0.357] | [3.6·10⁻⁶, 3.6·10⁻⁴] | [0.019, 1.9] | T_ps [0.00646, 0.516] | T_ps [0.00381, 0.305] | [9.57·10⁻⁶, 9.57·10⁻⁴] |
| Landsat, K = 10, Gaussian | T_ps [0.00446, 0.357] | T_ps [0.00446, 0.357] | T_ps [0.00446, 0.357] | T_ps [0.00446, 0.357] | [1.6·10⁻⁵, 0.0016] | [0.0364, 3.64] | T_ps [0.00646, 0.516] | T_ps [0.00381, 0.305] | [7.3·10⁻⁵, 0.0073] |
| BlattHierarchy, K = 5, Cauchy | [0.00689, 0.689] | [0.00604, 0.604] | [0.0049, 0.49] | [0.00466, 0.466] | [6.5·10⁻⁴, 0.065] | [0.0249, 2.49] | [0.0128, 1.28] | [0.00559, 0.559] | [0.0016, 0.16] |

### Bridge

| Graph | Base | OR filter | No repair | Scale before repair | MAD bandwidth | No gauge | q = 5 | q = 40 | Compiler defaults |
|---|---|---|---|---|---|---|---|---|---|
| Iris, K = 5, Gaussian | 1: 1.7·10⁻⁷ | 1: 2.1·10⁻⁵ | none | 1: 3·10⁻⁸ | 1: 4.9·10⁻⁴¹ | 1: 6.1·10⁻⁷ | 1: 1.7·10⁻⁷ | 1: 1.7·10⁻⁷ | 1: 1.7·10⁻⁴⁰ |
| Iris, K = 10, Gaussian | 1: 5.7·10⁻⁶ | 1: 2·10⁻⁴ | none | 1: 4.5·10⁻⁶ | 1: 8.4·10⁻³¹ | 1: 4·10⁻⁵ | 1: 5.7·10⁻⁶ | 1: 5.7·10⁻⁶ | 1: 5.9·10⁻³⁰ |
| Iris, K = 10, Cauchy | 1: 0.0067 | 1: 0.0058 | none | 1: 0.0065 | 1: 4.8·10⁻⁴ | 1: 0.047 | 1: 0.0067 | 1: 0.0067 | 1: 0.0034 |

### iris-k5-gauss

| Arm | T_fs, T_ps (seeds with a cliff, if not 12) | T_fs window: reading | ARI to 50/100: largest, midpoint | Full window: reading; seeds `ColdCensored` | [T_fs, T_ps]: reading | Species at T* (≥ 9 seeds) |
|---|---|---|---|---|---|---|
| Base | 0.00442, 0.038 | setosa + pair (12) | 0.965, 0.855 | setosa + pair (12); 12 | setosa + pair (12) | 0.0493–0.0638 |
| OR filter | 0.0611, 0.0826 | setosa + species (11), setosa + pair (1) | 0.293, 0.257 | setosa + pair (12); 12 | setosa + pair (12) | 0.0826 |
| No repair | 0.0515, 0.0826 | pieces + species (12) | 0.107, 0.084 | pieces + pair (12); 12 | pieces + species (12) | 0.0493–0.0826, setosa in pieces |
| Scale before repair | 0.00341, 0.038 | setosa + pair (12) | 0.965, 0.855 | setosa + pair (12); 12 | setosa + pair (12) | 0.0493–0.0638 |
| MAD bandwidth | 4.92·10⁻⁴, 0.00171 | setosa + species (11), setosa + one species (1) | 0.319, 0.255 | setosa + species (11), setosa + one species (1); 12 | setosa + species (12) | 4.71·10⁻⁴–0.00222 |
| No gauge | 0.0154, 0.127 | setosa + pair (12) | 0.965, 0.855 | setosa + pair (12); 12 | setosa + pair (12) | 0.179–0.232 |
| q = 5 | 0.00524, 0.0544 | setosa + pair (12) | 1.000, 0.938 | setosa + pair (12); 12 | setosa + pair (12) | 0.112–0.147 |
| q = 40 | 0.00411, 0.0325 | setosa + pair (12) | 0.965, 0.855 | setosa + pair (12); 12 | setosa + pair (12) | 0.0437 |
| Compiler defaults | 1.29·10⁻⁴, 0.00807 | setosa + pair (12) | 0.734, 0.683 | setosa + pair (12); 12 | setosa + pair (12) | 0.00287–0.00807 |
| No gauge, × K̂ lens | 0.0155, 0.133 | setosa + pair (12) | 0.965, 0.855 | setosa + pair (12); 12 | setosa + pair (12) | 0.172–0.223 |
| MAD bandwidth, `auto` | 4.59·10⁻⁴, 0.00164 | setosa + species (9), setosa + one species (3) | 0.319, 0.251 | setosa + species (9), setosa + one species (3); 12 | setosa + species (12) | 5.06·10⁻⁴–0.0022 |
| No gauge, `auto` | 0.0234, 0.143 | setosa + pair (12) | 0.947, 0.855 | setosa + pair (12); 12 | setosa + pair (12) | 0.313, setosa in pieces |
| Compiler defaults, `auto` | 1.54·10⁻⁴, 0.00813 | setosa + pair (8), setosa + species (4) | 0.720, 0.676 | setosa + pair (12); 12 | setosa + pair (8), setosa + species (4) | 0.00305–0.00737 |

### iris-k10-gauss

| Arm | T_fs, T_ps (seeds with a cliff, if not 12) | T_fs window: reading | ARI to 50/100: largest, midpoint | Full window: reading; seeds `ColdCensored` | [T_fs, T_ps]: reading | Species at T* (≥ 9 seeds) |
|---|---|---|---|---|---|---|
| Base | 0.00373, 0.0585 | setosa + pair (12) | 1.000, 0.815 | setosa + pair (12); 12 | setosa + pair (12) | 0.107 |
| OR filter | 0.0168, 0.09 | setosa + pair (12) | 1.000, 0.947 | setosa + pair (12); 12 | setosa + pair (12) | None |
| No repair | 0.0124, 0.0638 | setosa + pair (12) | 0.831, 0.791 | setosa + pair (12); 12 | setosa + pair (12) | 0.107 |
| Scale before repair | 0.00372, 0.0585 | setosa + pair (12) | 1.000, 0.815 | setosa + pair (12); 12 | setosa + pair (12) | None |
| MAD bandwidth | 10⁻⁴, 0.00442 | setosa + pair (12) | 0.748, 0.676 | setosa + pair (12); 12 | setosa + pair (12) | 0.00372–0.00482 |
| No gauge | none (0) | setosa + pair (12) | 1.000, 1.000 | setosa + pair (12); 12 | setosa + pair (12) | None |
| q = 5 | 0.00438, 0.0934 | setosa + pair (12) | 1.000, 0.839 | setosa + pair (12); 12 | setosa + pair (12) | None |
| q = 40 | 0.00377, 0.0475 | setosa + pair (12) | 1.000, 0.815 | setosa + pair (12); 12 | setosa + pair (12) | None |
| Compiler defaults | 0.00121, 0.0349 | setosa + pair (12) | 0.734, 0.668 | setosa + pair (12); 12 | setosa + pair (12) | 0.0294–0.038 |
| No gauge, × K̂ lens | 0.0285, 0.412 | setosa + pair (12) | 1.000, 0.815 | setosa + pair (12); 12 | setosa + pair (12) | 0.753 |
| MAD bandwidth, `auto` | 1.1·10⁻⁴, 0.00454 | setosa + pair (12) | 0.720, 0.676 | setosa + pair (12); 12 | setosa + pair (12) | 0.00373–0.00501 |
| No gauge, `auto` | 0.0442, 0.421 | setosa + pair (12) | 0.982, 0.815 | setosa + pair (12); 12 | setosa + pair (12) | None |
| Compiler defaults, `auto` | none (0) | setosa + pair (12) | 0.855, 0.734 | setosa + pair (12); 12 | setosa + pair (12) | None |

### iris-k10-cauchy

| Arm | T_fs, T_ps (seeds with a cliff, if not 12) | T_fs window: reading | ARI to 50/100: largest, midpoint | Full window: reading; seeds `ColdCensored` | [T_fs, T_ps]: reading | Species at T* (≥ 9 seeds) |
|---|---|---|---|---|---|---|
| Base | 0.0114, 0.0537 | setosa + pair (12) | 0.930, 0.783 | setosa + pair (12); 0 | setosa + pair (12) | 0.0826 |
| OR filter | 0.032, 0.0757 | setosa + pair (12) | 1.000, 0.896 | setosa + pair (12); 0 | setosa + pair (12) | None |
| No repair | 0.0092, 0.0585 | setosa + pair (12) | 0.896, 0.791 | setosa + pair (12); 12 | setosa + pair (12) | 0.0826 |
| Scale before repair | 0.0114, 0.0537 | setosa + pair (12) | 0.930, 0.783 | setosa + pair (12); 0 | setosa + pair (12) | None |
| MAD bandwidth | 0.00102, 0.0068 | setosa + pair (12) | 0.896, 0.783 | setosa + pair (12); 0 | setosa + pair (12) | 0.0135 |
| No gauge | none (0) | one lineage of 150 (12) | 0.000, 0.000 | one lineage of 150 (12); 12 | one lineage of 150 (12) | None |
| q = 5 | 0.0232, 0.0934 | setosa + pair (12) | 0.965, 0.783 | one lineage of 150 (12); 12 | setosa + pair (12) | None |
| q = 40 | 0.00878, 0.0437 | setosa + pair (12) | 0.965, 0.783 | setosa + pair (12); 0 | setosa + pair (12) | None |
| Compiler defaults | 0.0227, 0.107 | setosa + pair (12) | 1.000, 0.904 | setosa + pair (12); 0 | setosa + pair (12) | None |
| No gauge, × K̂ lens | 0.0769, 0.378 | setosa + pair (12) | 0.965, 0.783 | setosa + pair (12); 0 | setosa + pair (12) | 0.581 |
| Base, `auto` | 0.0109, 0.0522 | setosa + pair (12) | 0.965, 0.787 | setosa + pair (12); 12 | setosa + pair (12) | None |
| OR filter, `auto` | 0.0307, 0.0741 | setosa + pair (12) | 1.000, 0.909 | setosa + pair (12); 12 | setosa + pair (12) | None |
| No repair, `auto` | 0.0101, 0.0589 | setosa + pair (12) | 0.879, 0.799 | setosa + pair (12); 12 | setosa + pair (12) | None |
| Scale before repair, `auto` | 0.0117, 0.0511 | setosa + pair (12) | 0.960, 0.783 | setosa + pair (12); 12 | setosa + pair (12) | None |
| MAD bandwidth, `auto` | 0.00111, 0.00677 | setosa + pair (12) | 0.896, 0.800 | setosa + pair (12); 0 | setosa + pair (12) | None |
| No gauge, `auto` | 0.0766, 0.367 | setosa + pair (12) | 0.951, 0.783 | setosa + pair (12); 12 | setosa + pair (12) | 0.6 |
| q = 5, `auto` | 0.0223, 0.0881 | setosa + pair (12) | 0.982, 0.791 | setosa + pair (12); 12 | setosa + pair (12) | None |
| q = 40, `auto` | 0.00884, 0.0424 | setosa + pair (12) | 0.960, 0.783 | setosa + pair (12); 12 | setosa + pair (12) | None |
| Compiler defaults, `auto` | none (0) | setosa + pair (12) | 1.000, 1.000 | setosa + pair (12); 0 | setosa + pair (12) | None |

### landsat-k5-gauss

| Arm | T_fs, T_ps | Largest, 0.10: passing, test; purity; covered; classes | Midpoint, 0.10: same | Largest, passing at size rule, 0.10, 0.15, 0.25 | Midpoint, same | Full window: passing | [T_fs, T_ps], largest: passing; purity |
|---|---|---|---|---|---|---|---|
| Base | 0.0731, 0.0881 | 11, 12; 0.91; 2,391; 6 | 12, 12; 0.92; 1,720; 6 | 12, 11, 10, 9 | 12, 12, 12, 12 | 0 | 5; 0.80 |
| OR filter | 0.0881, 0.0967 | 0, 0; 0.59; 4,148; 4 | 0, 0; 0.58; 3,389; 4 | 0, 0, 0, 0 | 0, 0, 0, 0 | 0 | 0; 0.59 |
| No repair | 0.0881, 0.0967 | 12, 11; 0.90; 2,326; 6 | 12, 12; 0.91; 1,678; 6 | 12, 12, 7, 5 | 12, 12, 12, 12 | 0 | 3; 0.81 |
| Scale before repair | 0.0731, 0.0881 | 12, 12; 0.90; 2,520; 6 | 12, 11; 0.91; 1,721; 6 | 12, 12, 11, 7 | 12, 12, 12, 12 | 0 | 8; 0.86 |
| MAD bandwidth | 0.00446, 0.00537 | 12, 0; 0.93; 1,394; 5 | 12, 0; 0.95; 945; 5 | 12, 12, 12, 11 | 12, 12, 12, 12 | 12 | 11; 0.93 |
| No gauge | 0.224, 0.246 | 12, 12; 0.89; 2,646; 6 | 12, 12; 0.92; 2,036; 6 | 12, 12, 8, 8 | 12, 12, 12, 12 | 0 | 6; 0.84 |
| q = 5 | 0.106, 0.128 | 3, 3; 0.32; 5,428; 2 | 3, 3; 0.32; 4,924; 2 | 6, 3, 5, 5 | 7, 3, 5, 5 | 0 | 0; 0.32 |
| q = 40 | 0.0624, 0.0752 | 12, 11; 0.91; 2,314; 6 | 12, 12; 0.93; 1,602; 6 | 12, 12, 12, 9 | 12, 12, 12, 12 | 0 | 7; 0.86 |
| Compiler defaults | 0.0173, 0.0218 | 12, 0; 0.93; 1,661; 5 | 12, 0; 0.94; 1,162; 5 | 12, 12, 12, 12 | 12, 12, 12, 12 | 12 | 12; 0.93 |
| No gauge, × K̂ lens | 0.216, 0.248 | 12, 12; 0.90; 2,504; 6 | 12, 12; 0.92; 1,804; 6 | 12, 12, 9, 7 | 12, 12, 12, 12 | 0 | 6; 0.85 |
| MAD bandwidth, `auto` | 1.22·10⁻⁴, 3.6·10⁻⁴ (8) | 0, 0; 0.56; 4,501; 3 | 0, 0; 0.56; 4,266; 3 | 0, 0, 0, 0 | 0, 0, 0, 0 | 0 | 0; 0.56 |
| No gauge, `auto` | 0.221, 0.243 | 12, 12; 0.91; 2,576; 6 | 12, 12; 0.92; 1,958; 6 | 12, 12, 12, 7 | 12, 12, 12, 12 | 0 | 8; 0.87 |
| Compiler defaults, `auto` | none (0) | 0, 0; 0.32; 5,422; 2 | 0, 0; 0.32; 5,209; 2 | 0, 0, 0, 0 | 0, 0, 0, 0 | 0 | 0; 0.32 |

### landsat-k10-gauss

| Arm | T_fs, T_ps | Largest, 0.10: passing, test; purity; covered; classes | Midpoint, 0.10: same | Largest, passing at size rule, 0.10, 0.15, 0.25 | Midpoint, same | Full window: passing | [T_fs, T_ps], largest: passing; purity |
|---|---|---|---|---|---|---|---|
| Base | 0.0967, 0.117 | 9, 7; 0.94; 2,289; 6 | 2, 2; 0.64; 2,516; 4 | 4, 9, 9, 11 | 4, 2, 2, 10 | 0 | 0; 0.59 |
| OR filter | 0.0881, 0.106 | 1, 0; 0.58; 4,262; 4 | 1, 0; 0.59; 3,936; 4 | 1, 1, 1, 0 | 1, 1, 1, 1 | 0 | 0; 0.57 |
| No repair | 0.0967, 0.117 | 8, 6; 0.91; 2,438; 6 | 6, 5; 0.79; 2,300; 4 | 7, 8, 8, 9 | 7, 6, 6, 8 | 0 | 0; 0.59 |
| Scale before repair | 0.0967, 0.117 | 7, 6; 0.89; 2,202; 5 | 3, 2; 0.64; 2,394; 4 | 5, 7, 7, 9 | 5, 3, 3, 9 | 0 | 1; 0.59 |
| MAD bandwidth | 0.0094, 0.0113 | 12, 0; 0.93; 1,538; 5 | 12, 0; 0.94; 1,093; 5 | 12, 12, 12, 11 | 12, 12, 12, 12 | 12 | 12; 0.93 |
| No gauge | 0.27, 0.325 (7) | 0, 0; 0.25; 5,032; 1 | 0, 0; 0.33; 4,862; 2 | 0, 0, 0, 0 | 0, 0, 0, 0 | 0 | 0; 0.25 |
| q = 5 | 0.14, 0.169 | 9, 9; 0.93; 2,758; 6 | 2, 2; 0.66; 2,908; 4 | 11, 9, 9, 9 | 11, 2, 0, 8 | 0 | 0; 0.58 |
| q = 40 | 0.0752, 0.0906 | 7, 7; 0.88; 2,254; 5 | 4, 3; 0.62; 2,302; 4 | 7, 7, 7, 7 | 7, 4, 5, 6 | 0 | 0; 0.58 |
| Compiler defaults | 0.00446, 0.0094 | 0, 0; 0.55; 4,280; 3 | 0, 0; 0.58; 3,366; 4 | 0, 0, 0, 0 | 0, 0, 0, 0 | 0 | 0; 0.57 |
| No gauge, × K̂ lens | 0.514, 0.592 | 8, 6; 0.94; 2,113; 6 | 4, 4; 0.65; 2,358; 4 | 4, 8, 8, 10 | 4, 4, 4, 10 | 0 | 0; 0.59 |
| MAD bandwidth, `auto` | 6.62·10⁻⁴, 0.00145 | 0, 0; 0.57; 4,466; 3 | 0, 0; 0.57; 4,204; 3 | 0, 0, 0, 0 | 0, 0, 0, 0 | 0 | 0; 0.57 |
| No gauge, `auto` | 0.514, 0.625 | 9, 8; 0.94; 2,091; 6 | 3, 3; 0.60; 2,347; 4 | 3, 9, 9, 10 | 3, 3, 3, 8 | 0 | 0; 0.59 |
| Compiler defaults, `auto` | none (0) | 0, 0; 0.32; 5,722; 2 | 0, 0; 0.32; 5,410; 2 | 0, 0, 0, 0 | 0, 0, 0, 0 | 0 | 0; 0.32 |

### toy

| Arm | T_fs, T_ps | T_fs window, largest: test; clusters; ARI | Midpoint: same | Full window: test; clusters; ARI; seeds `ColdCensored` | [T_fs, T_ps]: test; clusters; ARI |
|---|---|---|---|---|---|
| Base | 0.0251, 0.0504 | 12; 3; 0.689 | 9; 3; 0.703 | 12; 4; 0.680; 12 | 12; 3; 0.692 |
| OR filter | 0.0607, 0.0881 | 12; 3; 0.688 | 12; 3; 0.695 | 0; 1; 0.000; 12 | 12; 3; 0.688 |
| No repair | 0.0251, 0.0504 | 12; 3; 0.691 | 10; 3; 0.704 | 12; 4; 0.680; 11 | 12; 3; 0.692 |
| Scale before repair | 0.0239, 0.0504 | 12; 3; 0.689 | 9; 3; 0.704 | 12; 4; 0.680; 12 | 12; 3; 0.690 |
| MAD bandwidth | 0.00446, 0.0059 | 0; 19; 0.088 | 0; 19; 0.066 | 0; 19; 0.088; 12 | 0; 18; 0.095 |
| No gauge | none (0) | 0; 1; 0.002 | 0; 1; 0.004 | 0; 1; 0.002; 12 | 0; 1; 0.002 |
| q = 5 | 0.0417, 0.0729 | 12; 3; 0.686 | 12; 3; 0.701 | 11; 4; 0.677; 6 | 12; 3; 0.686 |
| q = 40 | 0.0204, 0.0411 | 12; 3; 0.689 | 9; 3; 0.703 | 12; 4; 0.680; 11 | 12; 3; 0.694 |
| Compiler defaults | 0.00446, 0.00711 | 0; 20; 0.147 | 0; 20; 0.101 | 0; 20; 0.147; 0 | 0; 6; 0.493 |
| No gauge, × K̂ lens | 0.189, 0.398 | 12; 3; 0.692 | 11; 3; 0.704 | 12; 4; 0.680; 12 | 12; 3; 0.694 |
| MAD bandwidth, `auto` | 1.05·10⁻⁴, 2.79·10⁻⁴ | 0; 13; 0.409 | 0; 13; 0.386 | 0; 13; 0.409; 12 | 12; 3; 0.704 |
| No gauge, `auto` | 0.201, 0.4 | 12; 3; 0.689 | 7; 3; 0.704 | 12; 3; 0.679; 0 | 12; 3; 0.694 |
| Compiler defaults, `auto` | 3.58·10⁻⁴, 0.00209 | 2; 6; 0.586 | 0; 6; 0.583 | 2; 6; 0.586; 12 | 12; 3; 0.700 |

### blatt

| Arm | T_fs, T_ps | T_fs and full windows, largest: reading; clusters; ARI coarse/medium/fine | Midpoint ARI | Full: seeds `ColdCensored` | [T_fs, T_ps] | T* bands (≥ 11 seeds) |
|---|---|---|---|---|---|---|
| Base | 1.19·10⁻⁵, 2.08·10⁻⁵ | medium (12); 6; 0.332/1.000/0.353 | 0.332/1.000/0.353 | 0 | coarse (12); 2; 1.000/0.332/0.082 | coarse 5.36·10⁻⁶–2.87·10⁻⁵ ; medium 3.64·10⁻⁵–0.00347 ; fine 0.0044–0.0994 |
| OR filter | 1.1·10⁻⁵, 1.92·10⁻⁵ | medium (12); 6; 0.332/1.000/0.353 | 0.332/1.000/0.353 | 0 | coarse (12); 2; 1.000/0.332/0.082 | coarse 5.36·10⁻⁶–2.26·10⁻⁵ ; medium 2.87·10⁻⁵–0.00273 ; fine 0.0044–0.0994 |
| No repair | 10⁻⁶, 0.0277 | fine (12); 25; 0.066/0.291/0.884 | 0.066/0.291/0.884 | 12 | fine (12); 25; 0.066/0.291/0.884 | coarse – ; medium – ; fine 10⁻⁶–0.0782 |
| Scale before repair | 7.09·10⁻⁶, 1.19·10⁻⁵ | medium (12); 6; 0.332/1.000/0.353 | 0.332/1.000/0.353 | 0 | coarse (12); 2; 1.000/0.332/0.082 | coarse 3.32·10⁻⁶–1.4·10⁻⁵ ; medium 2.26·10⁻⁵–0.00215 ; fine 0.00273–0.0782 |
| MAD bandwidth | 1.38·10⁻⁶, 2.22·10⁻⁶ | medium (7), fine (5); 16; 0.205/0.729/0.533 | 0.200/0.717/0.515 | 0 | coarse (12); 2; 1.000/0.332/0.082 | coarse 10⁻⁶–2.61·10⁻⁶ ; medium 4.21·10⁻⁶–4.01·10⁻⁴ ; fine 5.09·10⁻⁴–0.03 |
| No gauge | 4.63·10⁻⁵, 7.48·10⁻⁵ | medium (12); 6; 0.332/1.000/0.353 | 0.332/1.000/0.353 | 0 | coarse (12); 2; 1.000/0.332/0.082 | coarse 2.26·10⁻⁵–9.51·10⁻⁵ ; medium 1.54·10⁻⁴–0.0115 ; fine 0.0146–0.329 |
| q = 5 | 1.78·10⁻⁵, 3.43·10⁻⁵ | medium (12); 6; 0.332/1.000/0.353 | 0.332/1.000/0.353 | 0 | coarse (12); 2; 1.000/0.332/0.082 | coarse 1.18·10⁻⁵–5.18·10⁻⁵ ; medium 8.48·10⁻⁵–0.00718 ; fine 0.00919–0.177 |
| q = 40 | 9.11·10⁻⁶, 1.71·10⁻⁵ | medium (12); 6; 0.332/1.000/0.353 | 0.332/1.000/0.353 | 0 | coarse (12); 2; 1.000/0.332/0.082 | coarse 5.24·10⁻⁶–2.17·10⁻⁵ ; medium 2.75·10⁻⁵–0.00247 ; fine 0.00397–0.068 |
| Compiler defaults | 4.39·10⁻⁶, 8.32·10⁻⁶ | fine (12); 21; 0.124/0.497/0.786 | 0.123/0.495/0.782 | 0 | coarse (12); 2; 1.000/0.332/0.082 | coarse 2.05·10⁻⁶–8.65·10⁻⁶ ; medium 1.4·10⁻⁵–0.00133 ; fine 0.00169–0.16 |
| No gauge, × K̂ lens | 4.48·10⁻⁵, 8.16·10⁻⁵ | medium (12); 6; 0.332/1.000/0.353 | 0.332/1.000/0.353 | 0 | coarse (12); 2; 1.000/0.332/0.082 | coarse 2.46·10⁻⁵–8.16·10⁻⁵ ; medium 1.32·10⁻⁴–0.0125 ; fine 0.0159–0.359 |
| Base, `auto` | 0.00689, 0.0246 | fine (12); 25; 0.080/0.344/0.985 | 0.071/0.313/0.928 | 12 | fine (12); 24; 0.082/0.352/0.997 | coarse – ; medium – ; fine 0.00689–0.0971 |
| OR filter, `auto` | 0.00604, 0.0575 | fine (12); 24; 0.082/0.353/1.000 | 0.082/0.352/0.998 | 12 | fine (12); 24; 0.082/0.353/1.000 | coarse – ; medium – ; fine 0.00604–0.0851 |
| No repair, `auto` | 0.0049, 0.0286 | fine (12); 26; 0.065/0.288/0.879 | 0.061/0.274/0.851 | 12 | fine (12); 25; 0.066/0.290/0.885 | coarse – ; medium – ; fine 0.0049–0.0691 |
| Scale before repair, `auto` | 0.00466, 0.0184 | fine (12); 25; 0.080/0.344/0.985 | 0.070/0.308/0.919 | 12 | fine (12); 24; 0.082/0.351/0.997 | coarse – ; medium – ; fine 0.00466–0.0656 |
| MAD bandwidth, `auto` | 6.5·10⁻⁴, 0.00256 | fine (12); 25; 0.080/0.345/0.985 | 0.069/0.305/0.913 | 12 | fine (12); 24; 0.082/0.353/0.998 | coarse – ; medium – ; fine 6.5·10⁻⁴–0.0297 |
| No gauge, `auto` | 0.0249, 0.089 | fine (12); 25; 0.080/0.344/0.985 | 0.071/0.312/0.926 | 12 | fine (12); 24; 0.082/0.352/0.998 | coarse – ; medium – ; fine 0.0249–0.351 |
| q = 5, `auto` | 0.0128, 0.0342 | fine (12); 25; 0.080/0.345/0.986 | 0.073/0.321/0.943 | 12 | fine (12); 24; 0.082/0.353/0.998 | coarse – ; medium – ; fine 0.0128–0.181 |
| q = 40, `auto` | 0.00559, 0.0221 | fine (12); 25; 0.080/0.344/0.984 | 0.071/0.311/0.924 | 12 | fine (12); 24; 0.082/0.352/0.997 | coarse – ; medium – ; fine 0.00559–0.0788 |
| Compiler defaults, `auto` | 0.0016, 0.00262 | fine (12); 24; 0.082/0.353/1.000 | 0.081/0.350/0.995 | 12 | fine (12); 21; 0.106/0.439/0.854 | coarse – ; medium – ; fine 0.00215–0.132 |

### theta

| θ | Iris K = 5: reading; ARI to 50/100 | Iris K = 10 Cauchy | Landsat K = 5, largest: passing, test; purity | Landsat K = 10 | Toy: test; clusters; ARI | BlattHierarchy: reading; ARI to medium |
|---|---|---|---|---|---|---|
| 0.3 | setosa + pair (12); 0.965 | setosa + pair (12); 1.000 | 9, 9; 0.90 | 5, 4; 0.65 | 12; 3; 0.685 | medium (12); 1.000 |
| 0.4 | setosa + pair (12); 0.965 | setosa + pair (12); 0.965 | 9, 11; 0.90 | 7, 7; 0.91 | 12; 3; 0.689 | medium (12); 1.000 |
| 0.5 | setosa + pair (12); 0.965 | setosa + pair (12); 0.930 | 11, 12; 0.91 | 9, 7; 0.94 | 12; 3; 0.689 | medium (12); 1.000 |
| 0.6 | setosa + pair (12); 0.965 | setosa + pair (12); 0.888 | 12, 10; 0.92 | 6, 6; 0.80 | 12; 3; 0.693 | medium (12); 1.000 |
| 0.7 | setosa + pair (12); 0.947 | setosa + pair (12); 0.888 | 12, 8; 0.91 | 7, 5; 0.91 | 10; 4; 0.619 | medium (12); 1.000 |

### species bands by theta

iris-k5-gauss-base 0.3 ['0.0638', '0.0826']
iris-k5-gauss-base 0.4 ['0.0638', '0.0826']
iris-k5-gauss-base 0.5 ['0.0493', '0.0638']
iris-k5-gauss-base 0.6 ['0.0493']
iris-k5-gauss-base 0.7 ['0.0493']
iris-k10-gauss-base 0.3 ['0.107']
iris-k10-gauss-base 0.4 ['0.107']
iris-k10-gauss-base 0.5 ['0.107']
iris-k10-gauss-base 0.6 []
iris-k10-gauss-base 0.7 []
