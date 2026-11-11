<!-- graphs -->
| Graph | Arm | Filter, K̂, repair edges, cut by the floor | Bandwidth rule, a | Gauge | J_first (nodes on its smaller side) | Lens: bottom, top, levels | χ stages | Shelves |
|---|---|---|---|---|---|---|---|---|
| Iris, K = 5, Gaussian | MAD bandwidth | mutual, 3.49, 15, 0 | MAD, 0.1212 | 0.2863 | 6.25·10⁻⁴ (38) | 5.13·10⁻⁵, 0.214, 91 | none | 0.0012, 0.00578 |
|  | No gauge | mutual, 3.49, 15, 0 | mean edge, 0.3065 | none | 6.06·10⁻⁷ (50) | 4.98·10⁻⁸, 2.25, 191 | 0.0173 | 1.89·10⁻⁵ |
|  | Compiler defaults | OR, 6.72, 1, 13 | MAD, 0.1212 | none | 0.00307 (38) | 2.52·10⁻⁴, 0.912, 89 | none | 0.00571, 0.0305 |
| Iris, K = 10, Gaussian | MAD bandwidth | mutual, 7.00, 5, 3 | MAD, 0.1414 | 0.1429 | 0.00203 (38) | 1.67·10⁻⁴, 0.205, 78 | none | 0.0128 |
|  | No gauge | mutual, 7.04, 5, 0 | mean edge, 0.3645 | none | 4.02·10⁻⁵ (50) | 3.3·10⁻⁶, 4.62, 153 | 0.0365 | 0.00141 |
|  | Compiler defaults | OR, 12.89, 1, 28 | MAD, 0.1414 | none | 0.0142 (38) | 0.00117, 1.54, 79 | none | 0.097 |
| Iris, K = 10, Cauchy | Base | mutual, 7.04, 5, 0 | mean edge, 0.3645 | 0.1421 | 0.00669 (50) | 5.49·10⁻⁴, 0.568, 76 | 0.0116 | 0.0128 |
|  | OR filter | OR, 13.27, 1, 0 | mean edge, 0.4754 | 0.0754 | 0.00584 (50) | 4.8·10⁻⁴, 0.48, 76 | 0.0303 | 0.0145 |
|  | No repair | mutual, 6.97, 0, 0 | mean edge, 0.3602 | 0.1434 | 0.0621 (38) | 0.0051, 0.57, 52 | none | none |
|  | MAD bandwidth | mutual, 7.04, 5, 0 | MAD, 0.09535 | 0.1421 | 4.78·10⁻⁴ (50) | 3.93·10⁻⁵, 0.12, 88 | 9.03·10⁻⁴ | 6.85·10⁻⁴ |
|  | Scale before repair | mutual, 7.04, 5, 0 | mean edge, 0.3602 | 0.1421 | 0.00653 (50) | 5.37·10⁻⁴, 0.565, 76 | 0.0115 | 0.0104 |
|  | No gauge | mutual, 7.04, 5, 0 | mean edge, 0.3645 | none | 0.0471 (50) | 0.00387, 4, 76 | 0.0714 | 0.0707 |
|  | q = 5 | mutual, 7.04, 5, 0 | mean edge, 0.3645 | 0.1421 | 0.00669 (50) | 9.33·10⁻⁴, 0.822, 74 | 0.022 | none |
|  | q = 40 | mutual, 7.04, 5, 0 | mean edge, 0.3645 | 0.1421 | 0.00669 (50) | 4.5·10⁻⁴, 0.485, 76 | 0.00771 | 0.00702 |
|  | Compiler defaults | OR, 13.27, 1, 0 | MAD, 0.09535 | none | 0.00337 (50) | 2.77·10⁻⁴, 1.05, 90 | 0.0196 | 0.0113, 0.165 |
| Toy, K = 10, Gaussian | MAD bandwidth | mutual, 7.55, 12, 581 | MAD, 0.05568 | 0.1325 | 3.47·10⁻⁷ (923) | 2.85·10⁻⁸, 0.188, 170 | 1.55·10⁻⁶ | 2.76·10⁻⁵ |
|  | No gauge | mutual, 7.91, 12, 0 | mean edge, 0.1622 | none | 0.00389 (53) | 3.19·10⁻⁴, 4.4, 104 | 0.189 | 0.228 |
|  | Compiler defaults | OR, 9.25, 0, 2053 | MAD, 0.05568 | none | 2.62·10⁻⁶ (930) | 2.15·10⁻⁷, 1.42, 170 | 1.41·10⁻⁵ | 1.73·10⁻⁴ |
| Landsat, K = 5, Gaussian | MAD bandwidth | mutual, 2.79, 983, 58 | MAD, 6.732 | 0.3579 | 2.05·10⁻¹⁰ (393) | 1.68·10⁻¹¹, 0.0408, 233 | 7.69·10⁻¹⁰, 1.21·10⁻⁴, 0.00274 | 1.67·10⁻⁴ |
|  | No gauge | mutual, 2.81, 983, 0 | mean edge, 22.71 | none | 0.154 (393) | 0.0127, 1.93, 55 | 0.207 | none |
|  | Compiler defaults | OR, 4.92, 0, 421 | MAD, 6.732 | none | 7.69·10⁻⁷ (436) | 6.32·10⁻⁸, 0.152, 159 | 7.09·10⁻⁷, 4.15·10⁻⁴, 0.0179 | 5.74·10⁻⁴, 0.026 |
| Landsat, K = 10, Gaussian | MAD bandwidth | mutual, 5.25, 282, 202 | MAD, 7.762 | 0.1903 | 4.8·10⁻⁶ (437) | 3.94·10⁻⁷, 0.0736, 132 | 4.38·10⁻⁶, 6.52·10⁻⁴, 0.00873 | 0.00165 |
|  | No gauge | mutual, 5.32, 282, 0 | mean edge, 23.84 | none | 0.326 (436) | 0.0267, 3.73, 54 | 0.527 | none |
|  | Compiler defaults | OR, 9.79, 0, 812 | MAD, 7.762 | none | 2.52·10⁻⁵ (440) | 2.07·10⁻⁶, 0.525, 135 | 2.54·10⁻⁵, 0.00366, 0.0747 | 0.00461 |
| BlattHierarchy, K = 5, Cauchy | Base | mutual, 3.62, 82, 0 | mean edge, 0.08313 | 0.2766 | 1.33·10⁻⁵ (600) | 1.1·10⁻⁶, 0.616, 144 | 1.11·10⁻⁵ | 9.22·10⁻⁶, 5.94·10⁻⁴, 0.0138 |
|  | OR filter | OR, 6.55, 20, 0 | mean edge, 0.1059 | 0.1526 | 1.2·10⁻⁵ (600) | 9.81·10⁻⁷, 0.51, 143 | 1.09·10⁻⁵ | 8.27·10⁻⁶, 4.88·10⁻⁴, 0.0218 |
|  | No repair | mutual, 3.48, 0, 0 | mean edge, 0.06391 | 0.2875 | 0.0388 (10) | 0.00319, 0.559, 57 | none | none |
|  | MAD bandwidth | mutual, 3.62, 82, 0 | MAD, 0.02741 | 0.2766 | 1.45·10⁻⁶ (600) | 1.19·10⁻⁷, 0.276, 159 | 1.33·10⁻⁶ | 1.06·10⁻⁶, 6.56·10⁻⁵, 0.00202 |
|  | Scale before repair | mutual, 3.62, 82, 0 | mean edge, 0.06391 | 0.2766 | 7.89·10⁻⁶ (600) | 6.48·10⁻⁷, 0.538, 148 | 6.9·10⁻⁶ | 6·10⁻⁶, 3.72·10⁻⁴, 0.0091 |
|  | No gauge | mutual, 3.62, 82, 0 | mean edge, 0.08313 | none | 4.83·10⁻⁵ (600) | 3.96·10⁻⁶, 2.23, 144 | 4.4·10⁻⁵ | 4.01·10⁻⁵, 0.00179, 0.05 |
|  | q = 5 | mutual, 3.62, 82, 0 | mean edge, 0.08313 | 0.2766 | 1.33·10⁻⁵ (600) | 1.86·10⁻⁶, 0.892, 142 | 1.81·10⁻⁵ | 7.74·10⁻⁴, 0.0239 |
|  | q = 40 | mutual, 3.62, 82, 0 | mean edge, 0.08313 | 0.2766 | 1.33·10⁻⁵ (600) | 8.99·10⁻⁷, 0.526, 144 | 1.1·10⁻⁵ | 8.75·10⁻⁶, 4.53·10⁻⁴, 0.0112 |
|  | Compiler defaults | OR, 6.55, 20, 0 | MAD, 0.02741 | none | 5.25·10⁻⁶ (600) | 4.31·10⁻⁷, 1.15, 160 | 4.84·10⁻⁶ | 3.34·10⁻⁶, 2.2·10⁻⁴, 0.0132 |

<!-- criteria -->
| Arm | 2. Seeds agree, of 63 (51) | 3. Landsat: passing on the four rows (9 each) | 4. Holds `auto` by the rule, of 63 (60) | Not a loss in the dataset's terms, of 63 | 5. Beats blind | 6. Longest list (4) |
|---|---|---|---|---|---|---|
| `chi up`, the primary | 52 | 12, 12, 12, 12 | 47 (strict 38, target 2, move 7) | 56 | Yes | 4 |
| `hot` | 52 | 12, 12, 12, 12 | 46 (strict 36, target 4, move 6) | 56 | – | 1 |
| `tfs` | – | 0, 0, 0, 0 | 43 (strict 37, target 0, move 6) | 56 | – | 1 |
| `chi band` | 52 | 12, 12, 12, 12 | 44 (strict 35, target 3, move 6) | 56 | Yes | 4 |
| `chi slice` | 52 | 4, 10, 12, 12 | 38 (strict 27, target 2, move 9) | 40 | Yes | 4 |
| `shelf up` | 41 | 0, 12, 0, 0 | 35 (strict 33, target 1, move 1) | 40 | Yes | 4 |
| `shelf band` | 41 | 0, 12, 0, 0 | 35 (strict 33, target 1, move 1) | 40 | Yes | 4 |
| `shelf slice` | 41 | 0, 12, 0, 0 | 37 (strict 35, target 1, move 1) | 40 | Yes | 4 |
| `null`, χ's length | – | 0, 12, 12, 12 | 40 (strict 38, target 1, move 1) | 44 | – | 4 |
| `null`, shelves' length | – | 0, 0, 0, 0 | 35 (strict 34, target 0, move 1) | 40 | – | 4 |

<!-- landsat -->
| Landsat graph | χ stages | Starts passing on 9 or more of 12, and their share of the lens | `auto` | `tfs` | `hot` |
|---|---|---|---|---|---|
| K = 5, MAD | 7.69·10⁻¹⁰, 1.21·10⁻⁴, 0.00274 | 7.43·10⁻⁴–0.00695, 9 of 78 | 0, 3 (3), 0.56, 70% | 0, 4 (3), 0.54, 65% | 12, 6 (5), 0.93, 28% |
| K = 5, compiler defaults | 7.09·10⁻⁷, 4.15·10⁻⁴, 0.0179 | 0.00445–0.0313, 8 of 53 | 0, 3 (2), 0.32, 84% | 0, 4 (3), 0.55, 67% | 12, 5 (5), 0.93, 25% |
| K = 10, MAD | 4.38·10⁻⁶, 6.52·10⁻⁴, 0.00873 | 0.00287–0.0152, 7 of 44 | 0, 3 (3), 0.57, 69% | 0, 4 (3), 0.54, 64% | 12, 5 (5), 0.93, 24% |
| K = 10, compiler defaults | 2.54·10⁻⁵, 0.00366, 0.0747 | 0.0203–0.108, 7 of 45 | 0, 2 (2), 0.32, 89% | 0, 4 (3), 0.55, 67% | 12, 5 (5), 0.94, 26% |

<!-- toy -->
| Toy graph | `auto` | `tfs` and `hot` | Full lens | At the shelf |
|---|---|---|---|---|
| MAD | 1 stripe on 12, ARI 0.41 | 3 stripes on 12, ARI 0.69 | 2 stripes on 12, ARI 0.68 | 2 stripes on 12, ARI 0.68 |
| No gauge | 3 stripes on 7, ARI 0.69 | 3 stripes on 8, ARI 0.69 | 0 stripes on 12, ARI 0.00 | 2 stripes on 12, ARI 0.68 |
| Compiler defaults | 2 stripes on 9, ARI 0.59 | 3 stripes on 12, ARI 0.69 | 1 stripe on 11, ARI 0.67 | 1 stripe on 10, ARI 0.68 |

<!-- blatt -->
| BlattHierarchy graph | `auto` | `hot` | Levels at ARI ≥ 0.8: `chi up` | `null`, shelves | `shelf band` | `shelf slice` |
|---|---|---|---|---|---|---|
| Base | fine on 12, passing 12, best ARI 0.98 | medium on 12, passing 12, best ARI 1.00 | –M– | –MF | CMF | CMF |
| OR filter | fine on 12, passing 12, best ARI 1.00 | medium on 12, passing 12, best ARI 1.00 | –M– | –MF | CMF | CMF |
| No repair | fine on 12, passing 12, best ARI 0.88 | fine on 12, passing 12, best ARI 0.88 | ––F | ––F | ––F | ––F |
| MAD bandwidth | fine on 12, passing 12, best ARI 0.99 | medium on 11, passing 3, best ARI 0.73 | ––– | ––F | CMF | CMF |
| Scale before repair | fine on 12, passing 12, best ARI 0.98 | medium on 12, passing 12, best ARI 1.00 | –M– | –MF | CMF | CMF |
| No gauge | fine on 12, passing 12, best ARI 0.98 | medium on 12, passing 12, best ARI 1.00 | –M– | –MF | CMF | CMF |
| q = 5 | fine on 12, passing 12, best ARI 0.99 | medium on 12, passing 12, best ARI 1.00 | –M– | –MF | –MF | –MF |
| q = 40 | fine on 12, passing 12, best ARI 0.98 | medium on 12, passing 12, best ARI 1.00 | –M– | –MF | CMF | CMF |
| Compiler defaults | fine on 12, passing 12, best ARI 1.00 | fine on 12, passing 3, best ARI 0.79 | ––– | ––F | C–F | CMF |
