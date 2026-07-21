# Modal Decomposition and Non-Intrusive Projection-Based Reduced-Order Model for Reservoir Pressure and Rate Transient Analysis

POD-based reduced-order modeling of reservoir pressure-transient response:
full-order simulation (MATLAB/MRST), modal decomposition and non-intrusive
ROM training/validation (Python), and inverse parameter characterization
against synthetic and published field data.

## Structure

```
case_studies/                  Ch.4 (FOM/PTA), Ch.5 (POD), Ch.6 (SSA)
                                five fixed reservoir cases, increasing complexity
  matlab/                       MRST single-case driver, one script per case
  notebooks/                    case_schematics, case1..case5, case2_mesh_comparison
  figures/                      ch4_*, ch5_*, ch6_*
  data/                         excluded from git, see data/README.md

ch07_non_intrusive_rom/        Ch.7 — parametric non-intrusive ROM (OpInf + RBF)
  matlab/                       batch training-data generation + single-case driver
  notebooks/                    ch07_rom_single_rate, ch07_rom_variable_rate
  rom/, rom_vrate/               trained operators (basis, interpolants)
  figures/                       ch7_*
  data/                         excluded from git, see data/README.md

ch08_inversion_characterization/  Ch.8 — inverse characterization (grid search,
                                   Nelder-Mead, two-stage, Horner cross-check)
                                   and external field-data validation
  matlab/                       validation-batch generation + single-case driver
  notebooks/                    ch08_inversion_characterization, ch08_examp_field_validation
  figures/                       ch8_*
  data/                         excluded from git, see data/README.md
```

Each notebook's first cell lists its own imports/dependencies.

## Data

Raw simulation data (~5.5GB: snapshot matrices, well data, rock/grid fields,
parameters) is hosted on Zenodo rather than in git:

- `case_studies/data/` — <ZENODO DOI LINK>
- `ch07_non_intrusive_rom/data/` — <ZENODO DOI LINK>
- `ch08_inversion_characterization/data/` — <ZENODO DOI LINK>

Extract into the matching `data/` folder so notebook-relative paths resolve.

## Environment

| Package | Version |
|---|---|
| MATLAB / MRST | R2023a / 2023a |
| Python | 3.11 |
| NumPy / SciPy | 1.26 / 1.12 |
| h5py / pandas | 3.10 / 2.1 |
