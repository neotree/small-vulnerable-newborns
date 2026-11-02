# Data and code availability

## Code

All analysis code is openly available in this repository:
`https://github.com/neotree/small-vulnerable-newborns`. The pipeline
reproduces every table, figure and number in the manuscript; see
`docs/MANUSCRIPT_MAP.md` for the script behind each result.

## Data

The individual-level neonatal and maternal data underlying this study were
collected during routine care through the Neotree platform at Kamuzu Central
Hospital (KCH, Malawi) and Sally Mugabe Central Hospital (SMCH, Zimbabwe).
They contain potentially identifiable clinical information and are held
under national and institutional data-governance agreements; they are not
publicly available and are not included in this repository. De-identified
data may be available from the corresponding author on reasonable request
and subject to the relevant data-governance and ethics approvals in Malawi
and Zimbabwe.

To reproduce the analysis with your own copy of the data, place the six
required files in `input/` exactly as described in `input/README.md`, then
run `Rscript run_all.R` from the repository root.

## Suggested statement for the manuscript

> Data and code availability. All analysis code is publicly available at
> https://github.com/neotree/small-vulnerable-newborns. The underlying
> individual-level data were collected during routine care via the Neotree
> platform and are governed by national and institutional data-sharing
> agreements in Malawi and Zimbabwe; they are not publicly available but may
> be requested from the corresponding author subject to the relevant
> approvals.
