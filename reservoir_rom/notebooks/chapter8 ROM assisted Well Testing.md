
======================================================================
 CASE 1/6: k=25.0 mD, S=5.0, sched=3
 Folder: k025_s05_sched3
======================================================================
  Source : F:\Niloy Deb\Part Three\reservoir_rom\notebooks\..\data\02_validation\k025_s05_sched3
  fwd_err: 0.0029  [EXCELLENT]

  Grid search 200x200=40000 pts ...
  Done 379.6s

  Method                   k [mD]   k%err       S   S err  calls
  --------------------------------------------------------
  True (FOM)                 25.0       —    5.00       —
  A: Grid Search            23.88    4.5%    4.42   0.578  40000
  C: Nelder-Mead           24.309   2.76%  4.6469  0.3531   3790
  B: Two-Stage              22.92     8.3%    3.90   1.101   ~110
  Horner (MTR)              23.55    5.8%   dp'=33.092 psi
  68% CI: k∈[23.9,23.9] S∈[4.42,4.42]  |  fwd_err: 0.0029 [EXCELLENT]

======================================================================
 CASE 2/6: k=45.0 mD, S=7.0, sched=9
 Folder: k045_s07_sched9
======================================================================
  Source : F:\Niloy Deb\Part Three\reservoir_rom\notebooks\..\data\02_validation\k045_s07_sched9
  fwd_err: 0.0014  [EXCELLENT]

  Grid search 200x200=40000 pts ...
  Done 420.5s

  Method                   k [mD]   k%err       S   S err  calls
  --------------------------------------------------------
  True (FOM)                 45.0       —    7.00       —
  A: Grid Search            45.19    0.4%    7.12   0.116  40000
  C: Nelder-Mead           44.538   1.03%  6.8745  0.1255    257
  B: Two-Stage              36.04    19.9%    3.98   3.020   ~110
  Horner (MTR)              42.18    6.3%   dp'=13.858 psi
  68% CI: k∈[45.2,45.2] S∈[7.12,7.12]  |  fwd_err: 0.0014 [EXCELLENT]

======================================================================
 CASE 3/6: k=85.0 mD, S=3.0, sched=4
 Folder: k085_s03_sched4
======================================================================
  Source : F:\Niloy Deb\Part Three\reservoir_rom\notebooks\..\data\02_validation\k085_s03_sched4
  fwd_err: 0.0026  [EXCELLENT]

  Grid search 200x200=40000 pts ...
  Done 411.2s

  Method                   k [mD]   k%err       S   S err  calls
  --------------------------------------------------------
  True (FOM)                 85.0       —    3.00       —
  A: Grid Search            83.78    1.4%    2.81   0.186  40000
  C: Nelder-Mead           85.850   1.00%  3.0816  0.0816   3809
  B: Two-Stage              92.85     9.2%    4.02   1.020   ~110
  Horner (MTR)              80.95    4.8%   dp'=8.024 psi
  68% CI: k∈[83.8,83.8] S∈[2.81,2.81]  |  fwd_err: 0.0026 [EXCELLENT]

======================================================================
 CASE 4/6: k=170.0 mD, S=5.0, sched=6
 Folder: k170_s05_sched6
======================================================================
  Source : F:\Niloy Deb\Part Three\reservoir_rom\notebooks\..\data\02_validation\k170_s05_sched6
  fwd_err: 0.0105  [EXCELLENT]

  Grid search 200x200=40000 pts ...
  Done 412.6s

  Method                   k [mD]   k%err       S   S err  calls
  --------------------------------------------------------
  True (FOM)                170.0       —    5.00       —
  A: Grid Search           202.93   19.4%    7.76   2.759  40000
  C: Nelder-Mead          205.178  20.69%  7.8377  2.8377    233
  B: Two-Stage             158.53     6.7%    3.98   1.020   ~110
  Horner (MTR)             159.15    6.4%   dp'=3.265 psi
  68% CI: k∈[202.9,202.9] S∈[7.76,7.76]  |  fwd_err: 0.0105 [EXCELLENT]

======================================================================
 CASE 5/6: k=220.0 mD, S=7.0, sched=9
 Folder: k220_s07_sched9
======================================================================
  Source : F:\Niloy Deb\Part Three\reservoir_rom\notebooks\..\data\02_validation\k220_s07_sched9
  fwd_err: 0.0125  [EXCELLENT]

  Grid search 200x200=40000 pts ...
  Done 206.1s

  Method                   k [mD]   k%err       S   S err  calls
  --------------------------------------------------------
  True (FOM)                220.0       —    7.00       —
  A: Grid Search           224.92    2.2%    7.64   0.638  40000
  C: Nelder-Mead          226.720   3.05%  7.7270  0.7270    107
  B: Two-Stage             179.36    18.5%    4.02   2.980   ~110
  Horner (MTR)             206.49    6.1%   dp'=2.831 psi
  68% CI: k∈[224.9,224.9] S∈[7.64,7.64]  |  fwd_err: 0.0125 [EXCELLENT]

======================================================================
 CASE 6/6: k=280.0 mD, S=5.0, sched=6
 Folder: k280_s05_sched6
======================================================================
  Source : F:\Niloy Deb\Part Three\reservoir_rom\notebooks\..\data\02_validation\k280_s05_sched6
  fwd_err: 0.0189  [EXCELLENT]

  Grid search 200x200=40000 pts ...
  Done 122.2s

  Method                   k [mD]   k%err       S   S err  calls
  --------------------------------------------------------
  True (FOM)                280.0       —    5.00       —
  A: Grid Search           293.89    5.0%    5.75   0.749  40000
  C: Nelder-Mead          292.365   4.42%  5.6738  0.6738    164
  B: Two-Stage             265.16     5.3%    4.10   0.899   ~110
  Horner (MTR)             264.89    5.4%   dp'=1.962 psi
  68% CI: k∈[293.9,293.9] S∈[5.75,5.75]  |  fwd_err: 0.0189 [EXCELLENT]

Batch done: 6/6 processed.
