# M9 relaxed torsional strain (2026-09-29T07:59:18)

tau_fold = max(|C3-C4-S1-N2|, |C5-C4-S1-N2|); strain = E(C3-C4-S1-N2 held at tau, everything else relaxed) - E(unconstrained PCM minimum), electronic, kcal/mol.

- sulfanilamide_pcm_d3: 22/22 points converged, reference: ref.json (unconstrained re-optimisation, converged=True); interpolation check (grid spline vs exact off-grid points) max |diff| 0.023 kcal/mol; mirror branch: 150 deg 0.005 kcal/mol
- sulfonamidate_pcm_d3: 18/18 points converged, reference: ref.json (unconstrained re-optimisation, converged=True); interpolation check (grid spline vs exact off-grid points) max |diff| 0.012 kcal/mol

| Species | Torsion | tau_fold (deg) | Relaxed strain (kcal/mol) | How | MMFF94s relaxed | Rigid DFT, burst 2 (superseded) |
|---|---|---:|---:|---|---:|---:|
| sulfanilamide | grid 90 | 90.0 | 0.01 | exact |  |  |
| sulfanilamide | grid 100 | 100.0 | 0.08 | exact |  |  |
| sulfanilamide | SAPH crystal (cation, CIF) | 102.7 | 0.13 | interpolated |  |  |
| sulfanilamide | 3TZF Y. pestis DHPS, bound sulfamethoxazole | 105.6 | 0.18 | exact |  |  |
| sulfanilamide | 3TZF (DHPS), Vina, representative (seed 2, Table S27) | 107.7 | 0.23 | exact | 0.10 | 3.33 |
| sulfanilamide | 3TZF (DHPS), Vinardo, representative (seed 1, Table S27) | 108.1 | 0.24 | interpolated | 0.10 | 3.35 |
| sulfanilamide | grid 110 | 110.0 | 0.28 | exact |  |  |
| sulfanilamide | 3TZF (DHPS), Vina, mean of the top poses | 110.5 | 0.28 | interpolated |  |  |
| sulfanilamide | 1AJ0 E. coli DHPS, bound sulfanilamide | 111.3 | 0.28 | exact |  |  |
| sulfanilamide | 7S2K (Sul2), Vina, representative (seed 4, Table S27) | 119.4 | 0.48 | interpolated | 0.28 | 3.56 |
| sulfanilamide | grid 120 | 120.0 | 0.50 | exact |  |  |
| sulfanilamide | 3TZF (DHPS), Vinardo, mean of the top poses | 124.6 | 0.61 | exact |  |  |
| sulfanilamide | 3HS4 (hCA II), Vinardo, mean of the top poses | 129.2 | 0.75 | exact |  |  |
| sulfanilamide | grid 130 | 130.0 | 0.77 | exact |  |  |
| sulfanilamide | 3HS4 (hCA II), Vinardo, representative (seed 1, Table S27) | 130.3 | 0.78 | interpolated | 0.47 | 3.88 |
| sulfanilamide | 6RL9 hCA II, bound sulfanilamide | 132.1 | 0.84 | exact |  |  |
| sulfanilamide | grid 140 | 140.0 | 1.08 | exact |  |  |
| sulfanilamide | 3HS4 (hCA II), Vina, representative (seed 2, Table S27) | 140.2 | 1.08 | interpolated | 0.85 | 4.54 |
| sulfanilamide | 3HS4 (hCA II), AD4Zn, mean of the top poses | 144.1 | 1.20 | exact |  |  |
| sulfanilamide | 3HS4 (hCA II), AD4Zn, representative (seed 13, Table S27) | 145.9 | 1.25 | interpolated | 1.05 | 4.75 |
| sulfanilamide | 3HS4 (hCA II), Vina, mean of the top poses | 146.3 | 1.26 | exact |  |  |
| sulfanilamide | 6RL9 (hCA II), AD4Zn, mean of the top poses | 147.0 | 1.29 | exact |  |  |
| sulfanilamide | 6RL9 (hCA II), AD4Zn, representative (seed 3, Table S27) | 148.6 | 1.34 | interpolated | 1.14 | 4.82 |
| sulfanilamide | grid 150 | 150.0 | 1.38 | exact |  |  |
| sulfanilamide | 6RL9 (hCA II), Vinardo, mean of the top poses | 156.0 | 1.53 | exact |  |  |
| sulfanilamide | 6RL9 (hCA II), Vinardo, representative (seed 9, Table S27) | 156.7 | 1.54 | interpolated | 1.44 | 5.02 |
| sulfanilamide | grid 160 | 160.0 | 1.63 | exact |  |  |
| sulfanilamide | 6RL9 (hCA II), Vina, representative (seed 9, Table S27) | 160.3 | 1.64 | interpolated | 1.56 | 5.07 |
| sulfanilamide | 6RL9 (hCA II), Vina, mean of the top poses | 160.4 | 1.64 | interpolated |  |  |
| sulfanilamide | 7S2K (Sul2), Vina, mean of the top poses | 162.5 | 1.70 | exact |  |  |
| sulfanilamide | grid 170 | 170.0 | 1.81 | exact |  |  |
| sulfanilamide | 7S2K (Sul2), Vinardo, representative (seed 1, Table S27) | 170.1 | 1.81 | interpolated | 1.79 | 5.15 |
| sulfanilamide | 7S2K (Sul2), Vinardo, mean of the top poses | 170.4 | 1.81 | interpolated |  |  |
| sulfanilamide | grid 180 | 180.0 | 1.89 | exact |  |  |
| sulfonamidate | grid 90 | 90.0 | -0.01 | exact |  |  |
| sulfonamidate | grid 100 | 100.0 | -0.00 | exact |  |  |
| sulfonamidate | SAPH crystal (cation, CIF) | 102.7 | 0.00 | interpolated |  |  |
| sulfonamidate | 3TZF Y. pestis DHPS, bound sulfamethoxazole | 105.6 | 0.00 | exact |  |  |
| sulfonamidate | grid 110 | 110.0 | 0.00 | exact |  |  |
| sulfonamidate | 1AJ0 E. coli DHPS, bound sulfanilamide | 111.3 | 0.00 | exact |  |  |
| sulfonamidate | 3TZF (DHPS), Vinardo, representative (seed 1, Table S27) | 113.7 | 0.00 | interpolated | 0.20 | 4.64 |
| sulfonamidate | 3TZF (DHPS), Vinardo, mean of the top poses | 113.8 | 0.00 | exact |  |  |
| sulfonamidate | 3TZF (DHPS), Vina, mean of the top poses | 115.0 | 0.00 | exact |  |  |
| sulfonamidate | 7S2K (Sul2), Vina, representative (seed 5, Table S27) | 115.2 | 0.00 | interpolated | 0.22 | 4.76 |
| sulfonamidate | 7S2K (Sul2), Vinardo, mean of the top poses | 115.3 | 0.00 | interpolated |  |  |
| sulfonamidate | 7S2K (Sul2), Vinardo, representative (seed 2, Table S27) | 115.3 | 0.00 | interpolated | 0.22 | 4.77 |
| sulfonamidate | 7S2K (Sul2), Vina, mean of the top poses | 115.5 | 0.00 | interpolated |  |  |
| sulfonamidate | 3TZF (DHPS), Vina, representative (seed 10, Table S27) | 115.7 | 0.00 | interpolated | 0.23 | 4.72 |
| sulfonamidate | grid 120 | 120.0 | -0.01 | exact |  |  |
| sulfonamidate | 3HS4 (hCA II), Vinardo, mean of the top poses | 130.0 | -0.04 | exact |  |  |
| sulfonamidate | grid 130 | 130.0 | -0.04 | exact |  |  |
| sulfonamidate | 3HS4 (hCA II), Vinardo, representative (seed 3, Table S27) | 132.0 | -0.05 | interpolated | 0.44 | 4.69 |
| sulfonamidate | 6RL9 hCA II, bound sulfanilamide | 132.1 | -0.05 | exact |  |  |
| sulfonamidate | grid 140 | 140.0 | -0.07 | exact |  |  |
| sulfonamidate | 3HS4 (hCA II), Vina, mean of the top poses | 149.4 | -0.09 | exact |  |  |
| sulfonamidate | grid 150 | 150.0 | -0.09 | exact |  |  |
| sulfonamidate | 3HS4 (hCA II), Vina, representative (seed 7, Table S27) | 150.4 | -0.09 | interpolated | 0.87 | 5.97 |
| sulfonamidate | 6RL9 (hCA II), Vinardo, mean of the top poses | 153.2 | -0.08 | exact |  |  |
| sulfonamidate | 6RL9 (hCA II), Vinardo, representative (seed 7, Table S27) | 153.6 | -0.08 | interpolated | 0.91 | 5.99 |
| sulfonamidate | 6RL9 (hCA II), Vina, mean of the top poses | 155.8 | -0.08 | exact |  |  |
| sulfonamidate | 6RL9 (hCA II), Vina, representative (seed 10, Table S27) | 155.8 | -0.08 | exact | 0.94 | 6.00 |
| sulfonamidate | grid 160 | 160.0 | -0.09 | exact |  |  |
| sulfonamidate | grid 170 | 170.0 | -0.09 | exact |  |  |
| sulfonamidate | grid 180 | 180.0 | -0.06 | exact |  |  |
