# Small and vulnerable newborns (SVN) analysis

R code for the analysis of small and vulnerable newborn phenotypes, stillbirth and neonatal mortality at Kamuzu Central Hospital (KCH, Malawi) and Sally Mugabe Central Hospital (SMCH, Zimbabwe), 1 January 2022 to 31 December 2025, using Neotree data. It produces every table, figure and number in the manuscript.

No patient data are included. The analysis runs on six locked Neotree exports held under the study's data governance arrangements (see `input/README.md`).

## Run

Use R 4.6.1 or later with the packages in `config/packages.csv` (`Rscript install_dependencies.R` installs any missing ones from CRAN). Run in a UTF-8 locale (for example `LC_ALL=C.UTF-8`). Place the six input files in `input/`, then:

```sh
Rscript run_all.R --check
Rscript run_all.R
```

`--check` verifies the inputs against `config/locked_inputs.csv`, the INTERGROWTH reference table and the packages. Each script runs first for all admissions (`output/all/`) and then for the ultrasound-dated subset (`output/uss/`); scripts with `uss_pass = FALSE` in `config/scripts.csv` run only once because they already compare the subgroups. To read inputs or write results elsewhere, set `LBW_INPUT_DIR` and `LBW_OUTPUT_DIR`.

`Rscript run_all.R --validate` runs the data-free checks used by GitHub Actions.

## Contents

- `scripts/`: 32 analysis scripts and the shared setup `00a_setup_shared.R`.
- `config/`: script order, input and reference checksums, packages and the antenatal-care helpers.
- `references/`: INTERGROWTH-21st 10th-centile birthweight table and its sources.
- `docs/MANUSCRIPT_MAP.md`: which script produces each manuscript item.
- `tests/`: data-free project checks.

## Classification rules

SGA is birthweight below the sex-specific INTERGROWTH-21st 10th centile at the recorded week and day (24+0 to 42+0); with unrecorded sex, the midpoint of the male and female centiles is used. Gestation `33.4` means 33 weeks 4 days. Invalid day endings are unclassified. Maternal births at exactly 22+0 and 23+0 use fetal estimated-weight 10th centiles (481 g, 538 g) as proxies. Preterm is below 37+0 weeks.
