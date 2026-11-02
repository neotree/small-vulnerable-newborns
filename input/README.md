# Inputs

No patient data are included. Place these six locked CSV exports here, or set `LBW_INPUT_DIR` to their folder. `run_all.R` checks each one against the MD5 in `config/locked_inputs.csv` and stops if any differs.

- `mwi-master-data-20220101-to-20251231.csv`
- `mwi-matched-nnu-data-20220101-to-20251231.csv`
- `mwi-maternal-data-20220101-to-20251231.csv`
- `zim-master-data-20220101-to-20251231.csv`
- `zim-matched-nnu-data-20220101-to-20251231.csv`
- `zim-maternal-data-20220101-to-20251231.csv`

`mwi` is KCH and `zim` is SMCH. Matched files hold neonatal admissions linked to discharge outcomes; master files hold every admission, including unmatched ones; maternal files hold delivery records. The Malawi files include `antenatalcare_category` (none, 1-3, >=4 visits), which is also filled where only a range of visits was recorded.
