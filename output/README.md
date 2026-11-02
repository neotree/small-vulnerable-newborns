# Outputs

`run_all.R` writes:

- `all/tables/` and `all/graphs/`: results for all admissions, including the full versus ultrasound versus non-ultrasound comparisons (S8, S9) and S10 (`all/tables/supplementary_s10/`).
- `uss/tables/` and `uss/graphs/`: the same analyses restricted to ultrasound-dated admissions.
- `logs/`: per-script logs (ultrasound pass in `logs/uss/`), input checksums, run summary and R session information.

`all/tables/manuscript_text_values.csv` lists every derived number quoted in the manuscript text with its source table. Re-running overwrites these files; set `LBW_OUTPUT_DIR` to keep runs separate.
