# Manuscript to R code map

Each manuscript item and the R script that produces it. Script names keep their internal numbering, which need not match the manuscript numbering.

| Manuscript item | R script |
|---|---|
| Table 1 and Results | `helper_Table01_SVN_by_site_and_stillbirth_OR.R` |
| Table 1 | `Table_01a_maternal_characteristics.R` |
| Table 1 | `Table_01b_stillbirth_rate_and_OR.R` |
| Supplementary Figure S1 and Results | `Figure_01_cohort_flowchart.R` |
| Supplementary Figure S2 | `Figure_02_stillbirth_OR_by_SVN.R` |
| Table 2 | `Table_02_NNU_characteristics.R` |
| Table 3 and Supplementary Table S1 | `Table_03_04_distribution_and_CFR.R` |
| Supplementary Tables S2a and S2b | `Table_05_timing_of_death.R` |
| Supplementary Figure S3 | `Figure_05_timing_of_death_by_SVN.R` |
| Supplementary Table S3 | `Table_06_length_of_stay.R` |
| Results | `helper_Table03_06_kruskal_wallis_cfr_los.R` |
| Supplementary Table S4a | `Table_07a_causes_of_death.R` |
| Supplementary Tables S4b and S4c | `Table_07bc_causes_of_death_by_site.R` |
| Table 4 and Results | `Table_08_univariable_predictors_AllSVN.R` |
| Supplementary Table S6 | `Table_09_temperature_categorical_ORs.R` |
| Figure 1 dependency | `dep_Table09_temperature_categorical_OR.R` |
| Results temperature interaction | `helper_Table09_temperature_svn_interaction.R` |
| Supplementary Table S5 | `Table_10_cfr_by_resuscitation.R` |
| Supplementary Table S7 and Results | `Table_11_CFR_by_antenatal_steroids.R` |
| Results and Discussion | `MV_regression_SVN_population_continuous_temperature.R` |
| Results and Discussion | `MV_regression_whole_cohort_continuous_temperature.R` |
| Results | `helper_HIV_CFR_exposed_vs_unexposed.R` |
| Supplementary Table S9 | `Appendix_Table_Sx_univariable_full_vs_USS_vs_nonUSS.R` |
| Supplementary Table S8 | `Supplementary_Table_S1_validation_CFR_full_vs_USS.R` |
| Figure 1 dependency and Supplementary Table S6 | `Figure_07_temperature_categorical_ORs_by_SVN.R` |
| Supplementary Table S6 | `helper_Table09_temperature_AllSVN_column.R` |
| Supplementary Table S10 and Results | `Supplementary_Table_S10_outcome_exclusion_comparison.R` |
| Table 4 continuous temperature and ANC Results | `Manuscript_temperature_and_ANC.R` |
| Figure 1 | `Figure_1_manuscript_temperature.R` |
| Table 3 and Supplementary Table S1 | `Manuscript_table_presentations.R` |
| Results (HIV distribution) | `helper_HIV_SVN_distribution.R` |
| Abstract, Methods, Results, Discussion and supplementary notes (derived text values) | `Manuscript_text_values.R` |

`config/scripts.csv` sets the execution order. Every script uses `scripts/00a_setup_shared.R` and the INTERGROWTH reference table. `Manuscript_text_values.R` writes every derived number quoted in the text, with its source table, to `manuscript_text_values.csv`. Ultrasound-dated results for S8, S9 and the related text come from `output/all/`; `output/uss/` holds the same analyses restricted to ultrasound-dated admissions.
