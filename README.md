# Reservoir ROM — MSc Thesis Code

Code accompanying my MSc thesis on non-intrusive reduced-order modeling (POD-based)
of reservoir pressure-transient response, including full-order simulation (MATLAB/MRST),
modal decomposition and ROM training/validation (Python), and parameter inversion
against both synthetic and published field data.

## Repository structure

```
reservoir_rom/                     Main pipeline: FOM data generation, POD analysis,
                                    SSA diagnostics, non-intrusive ROM training (Ch.7)
                                    and inversion characterization (Ch.8)
  matlab/                          MRST simulation drivers (single-case + batch)
  notebooks/                       POD / SSA / ROM / inversion notebooks
  figures*/                        Generated figures (constant-rate, variable-rate,
                                    inversion, external/example validation, ...)
  rom/, rom_vrate/, romshift/       Trained ROM operators (basis, interpolants)
  results/                         Prediction/validation outputs
  data/                            Raw simulation data — NOT in git, see data/README.md

reservoir_rom_field_validation/     External field-data validation pipeline
                                    (renamed from "reservoir_rom - for general
                                    validation - ofd check"). Notebooks 00-05 run
                                    simulation QC -> POD basis -> operator inference
                                    -> operator interpolation -> ROM validation ->
                                    ROM prediction, cross-checked against Horner
                                    analysis on the published well-test data below.
  data/                            Raw data — NOT in git, see data/README.md

well_test_reference_data/          Published pressure drawdown/buildup example
                                    datasets (Horne, 1990) used for external field
                                    validation — matches Appendix E of the thesis.

exploratory/reservoir_rom_CO/       Wider (k,S) parameter sweep that does NOT match
                                    the thesis's documented training/validation
                                    grids and isn't referenced anywhere in the
                                    thesis text. Kept for provenance only — not
                                    part of the reported results. See its own
                                    README for details.
```

## Data

Raw simulation data (~8GB total: snapshot matrices, well data, rock/grid fields,
parameter metadata) is hosted separately on Zenodo rather than in this repository,
to keep the git history lightweight:

- `reservoir_rom/data/` — <ZENODO DOI LINK — TODO>
- `reservoir_rom_field_validation/data/` — <ZENODO DOI LINK — TODO>
- `exploratory/reservoir_rom_CO/data/` — <ZENODO DOI LINK — TODO>

Download and extract each archive into the corresponding `data/` folder so that
notebook-relative paths resolve correctly. See the `README.md` inside each `data/`
folder for details.

## Notes on provenance

This repo was consolidated from three incremental work exports. Where the same
file existed in more than one export with different content, the most recently
modified version was kept (verified file-by-file, not by export order — the
later export was not always the most recently edited one). An earlier,
DMD-based (Dynamic Mode Decomposition) line of analysis was excluded entirely,
as the thesis's modal decomposition method is POD, not DMD.

## Software environment

| Package | Version | Language | Purpose |
|---|---|---|---|
| MATLAB | R2023a | MATLAB | Simulation driver |
| MRST | 2023a | MATLAB | Reservoir simulation framework |
| Python | 3.11 | Python | Modal decomposition, ROM, inversion |
| NumPy | 1.26 | Python | Array ops, SVD, linear algebra |
| SciPy | 1.12 | Python | Optimization, `.mat` I/O, interpolation |
| h5py | 3.10 | Python | HDF5 snapshot access |
| pandas | 2.1 | Python | Data tabulation |
