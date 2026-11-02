# INTERGROWTH reference data

## Official source tables

- [Very preterm boys, 24–32 weeks, birthweight centiles](https://intergrowth21.com/sites/default/files/2025-09/grow_verypreterm-ct-boys_bw_table.pdf) (Villar et al., *Lancet* 2016).
- [Very preterm girls, 24–32 weeks, birthweight centiles](https://intergrowth21.com/sites/default/files/2025-09/grow_verypreterm-ct-girls_bw_table.pdf) (Villar et al., *Lancet* 2016).
- [Newborn boys, 33–42 weeks, birthweight centiles](https://intergrowth21.com/sites/default/files/2023-02/grow_newborn-ct-boys_bw_table.pdf) (Villar et al., *Lancet* 2014).
- [Newborn girls, 33–42 weeks, birthweight centiles](https://intergrowth21.com/sites/default/files/2023-02/grow_newborn-ct-girls_bw_table.pdf) (Villar et al., *Lancet* 2014).
- [Fetal estimated weight, exact weeks, pooled sex](https://intergrowth21.com/sites/default/files/2023-01/grow_efw_ct_table_values.pdf) (Stirnemann et al., *Ultrasound in Obstetrics & Gynecology* 2017).

`references/official_intergrowth_p10_birthweight_24_42.csv` is the transcribed 10th-centile column of the four birthweight tables.
All 266 entries in this CSV were cross-checked against the 10th-centile rows in the four linked Oxford PDFs on 22 October 2026; no mismatches were found. The midpoint and fetal EFW proxies follow the agreed manuscript study rules. They are not official combined-sex birthweight standards.


The lookup is included as a required public reference input. The runner verifies its checksum in `config/reference_lock.csv`. Rows after 42+0 remain in the published-table transcription but are outside this manuscript analysis boundary.
