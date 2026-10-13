# Data availability note

Last updated: 2026-10-12. Sources downloaded 2026-10-07 to 2026-10-12; gnomAD counts pulled through the API on 2026-10-11.

| Source | Version / table | Used for | Access | Local file |
|---|---|---|---|---|
| gnomAD | v4.1.2, GRCh38 (dataset gnomad_r4) | Allele counts and frequencies by genetic ancestry group | API (`scripts/03_pull_gnomad.py`) and website tables | `data/gnomad_variant_frequencies.csv`, `data/gnomad_variant_frequencies_combined.csv`, `data/gnomAD_v4_samplecounts.csv` |
| CPIC | Gene-drug list (319 pairs); per-gene allele definition, phenotype and frequency tables | Level A gene-drug pairs, allele definitions, phenotypes | TSV and Excel downloads | `data/raw/cpic-genes-drugs.tsv`, `data/raw/cpic/` |
| ClinPGx (formerly PharmGKB) | Annotation downloads | Variant annotations and study populations | Downloads page | `data/raw/clinpgx/` |
| ClinVar | Records for the five selected variants | Clinical significance cross-check | Website | `data/variant_table.csv` |
| Statistics Canada | 2021 census, table 98-10-0351-01 | Population-group counts for Canada | CSV download | `data/raw/statcanada/98100351-eng/` |
| Ensembl VEP | Ensembl release 116 (June 2026); GRCh38; archive jun2026.archive.ensembl.org | Transcript consequences, HGVS coding/protein annotations, SIFT and PolyPhen predictions for five selected variants | Manual web-tool run; downloaded 2026-10-12 | `data/vep_results.txt`; processed using `scripts/04_vep_annotation.py` into `data/vep_annotations.csv` | 

### VEP annotation and identity check
Annotations come from the Ensembl VEP web tool (release 116, June 2026, GRCh38), run on 2026-10-12 with the MANE, HGVS, SIFT and PolyPhen options on. The raw output is saved in `data/vep_results.txt`, and `scripts/04_vep_annotation.py` turns it into `data/vep_annotations.csv`: one row per variant, using the MANE Select transcript. Cells with no value are written as "not applicable".

Identity check (2026-10-12): the rsID, GRCh38 position and (where VEP gives them) HGVS cDNA and protein names were compared with the ClinVar record and the CPIC allele definition table. ABCG2, CYP3A5, CYP4F2 and NUDT15 agree on all of these. VKORC1's variant lies upstream of the transcript, so VEP gives no HGVS name; it was matched on rsID, position and ClinVar's c.-1639G>A record. For ABCG2, VKORC1, CYP3A5 and CYP4F2 the genomic allele is the complement of the transcript (HGVS c.) allele, because these genes lie on the reverse strand; NUDT15 reads the same on both.

Consequences: missense for ABCG2 (Q141K), CYP4F2 (V433M) and NUDT15 (R139C); intron_variant for CYP3A5 `*3` (CPIC describes it as a splicing defect); upstream_gene_variant for VKORC1. SIFT and PolyPhen are computational predictions of protein impact, not evidence of drug response: ABCG2 Q141K and NUDT15 R139C are predicted tolerated and benign even though CPIC classes them as decreased and no-function alleles. This check covers variant identity and nomenclature only, not every ClinVar submission or CPIC recommendation.

## Key findings
- **gnomAD:** 807,162 individuals (730,947 exomes; 76,215 genomes). European groups are 77% and Middle Eastern 0.38%. In the summary table, Amish are counted within Remaining and Finns within European.
- **gnomAD counts for the five variants:** 86 rows by variant, dataset and group, and 50 combined rows (5 variants × 10 groups). AF equals AC/AN in both files, the combined AN equals exome AN plus genome AN, and one row (ABCG2, Middle Eastern) was checked by hand.
- **Coverage:** the VKORC1 variant is non-coding and appears only in the genomes. Amish are genomes-only (AN about 910 for every variant, and NUDT15 has no alternate alleles observed). The lowest AN is VKORC1 in Middle Eastern (292).
- **Precision flags** (in `data/gnomad_variant_frequencies_combined.csv`): 7 rows `low_an` (five Amish rows, VKORC1 Middle Eastern, VKORC1 Remaining), 1 row `few_alt` and `not_observed` (NUDT15 Amish).
- **CPIC:** 93 Level A gene-drug pairs across 21 genes (levels: A 93, B 15, C 207, retired 4). This project uses Level A only. Allele definition, phenotype and frequency tables are per-gene Excel downloads from the CPIC gene pages.
- **Selected variants:** five variants, one per gene (ABCG2, VKORC1, CYP3A5, CYP4F2, NUDT15); see `data/variant_table.csv`.
- **Panel (Step 5):** the all-five panel is at or above 99% in eight of ten groups (lowest: African/African American 69%, South Asian 98.7%), so the reduced panel (without CYP3A5 and VKORC1) is the more informative summary; see `data/panel_by_group.csv`.
- **VEP:** the five variants were annotated on their MANE Select transcripts (three missense, one intronic, one upstream); see the VEP section above.
- **ClinPGx:** `variantAnnotations/study_parameters.tsv` has 36,083 rows; population is in `Biogeographical Groups`, and about 27% are Unknown. European and East Asian are the largest named groups.
- **Statistics Canada:** Canada total 36,328,475 (persons in private households, 25% sample data), in 13 categories. Indigenous people are not a separate category and fall under "Not a visible minority". Canada's row excludes some incompletely enumerated reserves.

## Decisions
1. **gnomAD data.** The full gnomAD v4 release is used (not the non-UKB subset). Primary: exomes and genomes combined (AC and AN added). Sensitivity checks: exomes only (four genes; the VKORC1 variant is not in the exome data) and genomes only (all five genes, so every gene uses the same dataset). VKORC1 therefore uses genomes in the primary analysis too. Exome and genome samples differ in composition and sequencing, so differences between them are read cautiously.
2. **What is counted.** For each gene, the share of people carrying at least one copy of the variant allele (gnomAD's alternate allele), calculated as 1 − (1 − AF)². For CYP3A5, a second estimate gives expressers (at least one working copy), 1 − AF(`*3`)², counting `*3` as the only non-working allele; other no-function alleles (`*6`, `*7`) are ignored, so expressers are somewhat overestimated, most in African ancestry groups. These are carrier estimates, not the share needing a prescribing change. The five-gene panel (at least one of five) is shown with and without the two most common variants (CYP3A5 `*3` and VKORC1 A allele), because it is at or above 99% in eight of the ten groups. Genotype-to-phenotype mappings and recommendations are to be checked against the current CPIC guidelines.
3. **gnomAD groups.** All 10 groups are kept with their original labels. Small samples and low or missing allele numbers are flagged, not excluded. Mapping to census groups is documented separately.
4. **Precision flags.** Rows (variant × ancestry group) with an allele number below 2,500 are flagged `low_an`: at this AN the approximate 95% margin of error is about ±2 percentage points even at a frequency of 50%, the worst case. Rows with fewer than 20 observed alternate alleles (AC < 20) are flagged `few_alt`: at that count the approximate 95% margin of error is about ±44% of the estimated frequency. AC = 0 means the allele was not observed, not that its frequency is zero. A missing or zero AN is stored as blank. Flagged rows are kept and marked in the results, and are dropped in a sensitivity check. Confidence intervals use the Wilson method.
5. **Annotation.** Variants are annotated on the MANE Select transcript.
6. **Panel estimates.** The probability of carrying at least one of the selected variants is 1 − ∏(1 − cᵢ), where cᵢ is each gene's carrier probability. There is one variant per gene, so nothing is combined within a gene, and NUDT15's three drugs share one variant and are counted once. The product assumes the variants are independent. The five variants lie on different chromosomes (4, 7, 13, 16 and 19), so they are not linked, but population structure and mixed ancestry can still correlate them within a group, so the combined probability is an approximation. Panel bounds combine the lower and upper gene-level bounds, which makes them conservative. The panel is reported two ways: all five genes, and without CYP3A5 and VKORC1, the two genes with the highest median carrier probability across groups (0.99 and 0.63, against 0.48, 0.19 and 0.02 for CYP4F2, ABCG2 and NUDT15). The all-five panel is at or above 99% in eight of the ten groups, so it separates groups poorly and the reduced panel is the more informative summary. Panel rows are flagged when any contributing gene row is flagged (3 of 10 groups for the all-five panel, 1 of 10 for the reduced panel). The panel describes carriage of at least one of the five variants, not the proportion of people who would need a prescribing change.

| Gene | gnomAD alternate allele | Same as |
|---|---|---|
| ABCG2 | T | Q141K |
| VKORC1 | T | -1639G>A (A allele on the gene's strand) |
| CYP3A5 | C | `*3` |
| CYP4F2 | T | V433M (`*3`) |
| NUDT15 | T | `*3` (R139C) |

## Mapping issues
gnomAD (genetic ancestry groups), the census (self-reported population groups), ClinPGx (biogeographical groups) and CPIC's frequency tables (their own group labels) don't match each other. A mapping table and a sensitivity analysis on it will be needed. No gnomAD group corresponds to Indigenous peoples.

## Still open
- Counting unit for the evidence count (rows, variant annotations or papers; studies or participants)
- Evidence base: ClinPGx-annotated literature or the papers cited in CPIC guidelines
- Handling of Indigenous peoples in the census table (exclude, or separate them using table 98-10-0324-01)

## Census-to-gnomAD mapping

| Data file | Generated by | Last updated |
| --- | --- | --- |
| `data/census_gnomad_mapping.csv` | `scripts/07_census_mapping.py` | 2026-10-13 |
