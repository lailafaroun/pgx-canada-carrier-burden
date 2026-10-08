# Data availability note

Last updated: 08-Oct-2026

| Source | Version / table | Used for | Access | Local file |
|---|---|---|---|---|
| gnomAD | v4.1.2, GRCh38 | Allele frequencies by genetic ancestry group | Website tables | `data/gnomAD_v4_samplecounts.csv` |
| CPIC | Gene-drug list, 319 pairs | Level A gene-drug pairs; allele definition tables (per-gene Excel files, downloaded at variant selection) | TSV and Excel downloads | `data/raw/cpic-genes-drugs.tsv` |
| ClinPGx (formerly PharmGKB) | Annotation downloads | Variant annotations and study populations | Downloads page | `data/raw/clinpgx/` |
| Statistics Canada | 2021 census, table 98-10-0351-01 | Population-group counts for Canada | CSV download | `data/raw/statcanada/98100351-eng/` |

## Key findings
- **gnomAD:** 807,162 individuals (730,947 exomes; 76,215 genomes). European groups are 77% and Middle Eastern 0.38%. In the summary table, Amish are counted within Remaining and Finns within European.
- **CPIC:** 93 Level A gene-drug pairs across 21 genes (levels: A 93, B 15, C 207, retired 4). This project uses Level A only.
- **ClinPGx:** `variantAnnotations/study_parameters.tsv` has 36,083 rows; population is in `Biogeographical Groups`, and about 27% are Unknown. European and East Asian are the largest named groups.
- **Statistics Canada:** Canada total 36,328,475 (persons in private households, 25% sample data), in 13 categories. Indigenous people are not a separate category and fall under "Not a visible minority". Canada's row excludes some incompletely enumerated reserves.

## Mapping issues
gnomAD (genetic ancestry groups), the census (self-reported population groups) and ClinPGx (biogeographical groups) don't match each other. A mapping table and a sensitivity analysis on it will be needed. No gnomAD group corresponds to Indigenous peoples.

## Open decisions
- gnomAD: exome, genome or combined frequencies; full dataset or non-UKB subset
- Counting unit for the evidence count (rows, variant annotations or papers; studies or participants)
- Evidence base: ClinPGx-annotated literature or the papers cited in CPIC guidelines
- Handling of Indigenous peoples in the census table (exclude, or separate them using table 98-10-0324-01)