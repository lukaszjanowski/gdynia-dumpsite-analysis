# Gdynia dumpsite MBES analysis repository

Reproducible R workflow supporting the manuscript analysis of bathymetric change at the Gdynia dredged-material disposal site. The production modules are derived from the verified project scripts supplied by the authors. Unmodified source copies are retained in `archive/source_scripts/`.

## Important year convention

The final survey was acquired during August-December **2023**. All publication-facing names use 2023. The archive preserves historical technical names containing 2024.

## Quick start

1. Open `gdynia-dumpsite-analysis-v2.Rproj`.
2. Place private inputs in `data/raw/` using names specified in `config/config.yml`.
3. Install packages listed in `DESCRIPTION`.
4. Run `source("run_all.R", encoding="UTF-8")`.

## Core workflow

- `01_inventory.R`: validates inputs and creates `outputs/paths.rds`.
- `02_harmonise_dems.R`: invokes the verified combined harmonisation and sequential co-registration source.
- `03_vertical_coregistration.R`: validates co-registration diagnostics as a distinct pipeline checkpoint.
- `04_dod_and_volumes.R`: raw DoDs, LoD scenarios and zonal volumes.
- `05_morphometry_acoustics.R`: terrain metrics and interpretive acoustic alignment.
- `06_depositional_patches.R`: connected components and dominant fields.
- `07_boundary_accounting.R`: exact cell-based allocation.
- `08_sensitivity_analysis.R`: threshold and minimum-area scenarios.
- `09_publication_figures.R`: Figures 1-6 and supplementary graphics.
- `10_publication_tables.R`: manuscript tables.
- `11_reproducibility_report.R`: session and output inventories.

## Figure sources

- Figure 1: asymmetric regional context.
- Figure 2: illuminated bathymetric DEMs.
- Figure 3: common-scale DoD intervals.
- Figure 4: two dominant fields, legend below map.
- Figures 5-6 and S1-S2: quantitative and sensitivity products from polished v4.

## Data policy

Raw bathymetric and acoustic data are ignored by Git. Commit metadata, code, small tables and permitted derived products only.

## Note on stages 02 and 03

The supplied and verified algorithm performs grid harmonisation and sequential vertical co-registration in one stateful script. V2 preserves this tested execution as `R/02_03_harmonise_coregister.R`; stages 02 and 03 remain separate public checkpoints without duplicating the expensive calculation.
