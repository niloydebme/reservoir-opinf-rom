# Reservoir ROM — MSc Thesis Code

Code accompanying my MSc thesis on non-intrusive reduced-order modeling (POD-based)
of reservoir pressure-transient response, including full-order simulation (MATLAB/MRST),
modal decomposition and ROM training/validation (Python), and parameter inversion
against both synthetic and published field data.

## Repository structure

```
case_studies/                      Ch.4 (FOM/PTA), Ch.5 (POD), Ch.6 (SSA) — five
                                    fixed reservoir cases of increasing complexity
  matlab/                          MRST single-case drivers, one per case
  notebooks/                       case_schematics + case1..case5 (cleaned to
                                    POD-only; see note below), case2_mesh_comparison
  figures/                         ch4_*, ch5_*, ch6_* — matches thesis naming 1:1
  data/                            Raw simulation data — NOT in git, see data/README.md

reservoir_rom/                     Main pipeline: non-intrusive ROM training (Ch.7)
                                    and inversion characterization (Ch.8)
  matlab/                          MRST simulation drivers (single-case + batch)
  notebooks/                       POD / SSA / ROM / inversion notebooks
  figures*/                        Generated figures (constant-rate, variable-rate,
                                    inversion, external/example validation, ...)
  rom/, rom_vrate/, romshift/       Trained ROM operators (basis, interpolants)
  results/                         Prediction/validation outputs
  data/                            Raw simulation data — NOT in git, see data/README.md
  (still being curated at the notebook level — see Status below)

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

exploratory/                       Work confirmed NOT part of the reported thesis
                                    results, kept for provenance only:
  reservoir_rom_CO/                 wider (k,S) parameter sweep, grid doesn't match
                                    the thesis's documented training/validation grids
  unused_case_variants/             e.g. case5_nonstationary.m (rate change + skin
                                    damage + infill well) — not referenced in the thesis
```

## Status

This repo is being curated in batches, going beyond simple file deduplication to
also prune notebook cells down to what the thesis actually reports, fix figure
names to match the thesis 1:1, and reorganize directories by chapter. Done so far:
`case_studies/` (fully cleaned — see note below). Still in progress: splitting
`reservoir_rom/` into `rom_single_rate/`, `rom_variable_rate/`, and
`inversion_characterization/`, and a lighter cleanup pass on
`reservoir_rom_field_validation/`.

**Note on `case_studies/`:** these 7 notebooks were recovered from an earlier
export where they'd been misfiled as DMD (Dynamic Mode Decomposition) content
and nearly excluded. Inspection showed they're almost entirely POD analysis
(DMD is mentioned only in a couple of legacy title/comment lines, now removed)
— they are the actual source of the Ch.4/5/6 figures for the five fixed cases.
Each notebook has been pruned of dead-end/superseded plotting cells (verified
via which `savefig()` calls were live vs. commented out, and which target
filenames were later overwritten by a refined cell), and renamed to match the
thesis's chapter/figure numbering.

## Data

Raw simulation data (~8.5GB total: snapshot matrices, well data, rock/grid fields,
parameter metadata) is hosted separately on Zenodo rather than in this repository,
to keep the git history lightweight:

- `case_studies/data/` — <ZENODO DOI LINK — TODO>
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
later export was not always the most recently edited one). A line of analysis
initially thought to be DMD-based (Dynamic Mode Decomposition, not the thesis's
method) turned out to be misfiled POD work and was recovered — see `case_studies/`
above.

Several notebooks still hardcode an absolute local path (e.g.
`Path(r'F:\Niloy Deb\Part One\...')`) from the original working environment —
this needs updating to a relative path before the notebooks will run elsewhere;
tracked as a follow-up.

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
