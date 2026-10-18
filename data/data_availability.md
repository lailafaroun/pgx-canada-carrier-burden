# Data availability note

Last updated: 2026-10-18. Sources downloaded 2026-10-07 to 2026-10-12; gnomAD counts pulled through the API on 2026-10-11.

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

- **Census mapping (Step 6):** S0 maps 34,476,095 people (94.90% of the census total). S1 excludes “Not a visible minority” and covers 7,786,820 people (21.43%); it is not a national estimate. S2 maps Filipino and Southeast Asian groups to East Asian, increasing coverage to 98.61%.
- **National estimates (Step 7):** Under S0, the all-five-gene panel estimate is 98.4%, the reduced panel estimate is 61.7%, NUDT15 carrier prevalence is 2.7%, and CYP3A5 expresser prevalence is 21.5%. These are weighted estimates for the mapped population under the stated assumptions.
- **Sensitivity analysis (Step 8):** Six scenarios (S0–S5) are available, with eight quantities per scenario (48 rows total). S3 excludes flagged population groups per estimate and renormalizes the remaining weights. S4 uses exome-only frequencies and S5 genome-only frequencies; VKORC1 is unavailable in exomes.
- **Interpretation:** The S1-excluded ranges are more informative about variation across the primary mapping and other sensitivity checks than the full ranges, because S1 omits the largest census category. These are scenario ranges, not confidence intervals.
- **Limitations:** Intervals reflect group-level gnomAD sampling uncertainty only. They do not capture census uncertainty, ancestry-mapping uncertainty, or assumptions in the carrier and panel models. The selected census table does not separately identify Indigenous people within “Not a visible minority.”
- See `data/national_estimates.csv`, `data/all_scenario_estimates.csv`, and `data/sensitivity_summary.csv` for results.

## Decisions
1. **gnomAD data.** The full gnomAD v4 release is used (not the non-UKB subset). Primary: exomes and genomes combined (AC and AN added). Sensitivity checks: exomes only (four genes; the VKORC1 variant is not in the exome data) and genomes only (all five genes, so every gene uses the same dataset). VKORC1 therefore uses genomes in the primary analysis too. Exome and genome samples differ in composition and sequencing, so differences between them are read cautiously.
2. **What is counted.** For each gene, the share of people carrying at least one copy of the variant allele (gnomAD's alternate allele), calculated as 1 − (1 − AF)². For CYP3A5, a second estimate gives expressers (at least one working copy), 1 − AF(`*3`)², counting `*3` as the only non-working allele; other no-function alleles (`*6`, `*7`) are ignored, so expressers are somewhat overestimated, most in African ancestry groups. These are carrier estimates, not the share needing a prescribing change. The five-gene panel (at least one of five) is shown with and without the two most common variants (CYP3A5 `*3` and VKORC1 A allele), because it is at or above 99% in eight of the ten groups. Genotype-to-phenotype mappings and recommendations are to be checked against the current CPIC guidelines.
3. **gnomAD groups.** All 10 groups are kept with their original labels. Small samples and low or missing allele numbers are flagged, not excluded. Mapping to census groups is documented separately.
4. **Precision flags.** Rows (variant × ancestry group) with an allele number below 2,500 are flagged `low_an`: at this AN the approximate 95% margin of error is about ±2 percentage points even at a frequency of 50%, the worst case. Rows with fewer than 20 observed alternate alleles (AC < 20) are flagged `few_alt`: at that count the approximate 95% margin of error is about ±44% of the estimated frequency. AC = 0 means the allele was not observed, not that its frequency is zero. A missing or zero AN is stored as blank. Flagged rows are kept and marked in the results, and are dropped in a sensitivity check. Confidence intervals use the Wilson method.
5. **Annotation.** Variants are annotated on the MANE Select transcript.
6. **Panel estimates.** The probability of carrying at least one of the selected variants is 1 − ∏(1 − cᵢ), where cᵢ is each gene's carrier probability. There is one variant per gene, so nothing is combined within a gene, and NUDT15's three drugs share one variant and are counted once. The product assumes the variants are independent. The five variants lie on different chromosomes (4, 7, 13, 16 and 19), so they are not linked, but population structure and mixed ancestry can still correlate them within a group, so the combined probability is an approximation. Panel bounds combine the lower and upper gene-level bounds, which makes them conservative. The panel is reported two ways: all five genes, and without CYP3A5 and VKORC1, the two genes with the highest median carrier probability across groups (0.99 and 0.63, against 0.48, 0.19 and 0.02 for CYP4F2, ABCG2 and NUDT15). The all-five panel is at or above 99% in eight of the ten groups, so it separates groups poorly and the reduced panel is the more informative summary. Panel rows are flagged when any contributing gene row is flagged (3 of 10 groups for the all-five panel, 1 of 10 for the reduced panel). The panel describes carriage of at least one of the five variants, not the proportion of people who would need a prescribing change.
7. **Census mapping.** Census groups are mapped to gnomAD groups as described in the census mapping section below. The primary mapping assigns "Not a visible minority" to European (non-Finnish) and covers 94.90% of the census total; S1 excludes "Not a visible minority" (21.43% coverage) and is a sensitivity scenario, not a national estimate.
8. **National estimate.** National values are weighted averages of the group-level values, using the census weights of the mapped gnomAD groups (the weights sum to 1 within each scenario). For the panels, the group panel values from Step 5 are averaged; national gene values are not recombined. Intervals are the weighted group-level Wilson bounds, so they reflect gnomAD sampling uncertainty only, not census rounding or sampling, the group mapping, or the Hardy-Weinberg and independence assumptions. Each estimate reports the share of weight that falls on flagged rows. Results are in `data/national_estimates.csv`, from `scripts/08_national_estimate.py`.

| Gene | gnomAD alternate allele | Same as |
|---|---|---|
| ABCG2 | T | Q141K |
| VKORC1 | T | -1639G>A (A allele on the gene's strand) |
| CYP3A5 | C | `*3` |
| CYP4F2 | T | V433M (`*3`) |
| NUDT15 | T | `*3` (R139C) |

## Mapping issues
gnomAD (genetic ancestry groups), the census (self-reported population groups), ClinPGx (biogeographical groups) and CPIC's frequency tables (their own group labels) don't match each other. The mapping and its sensitivity analysis are described in the census mapping section below. No gnomAD group corresponds to Indigenous peoples.

## Census-to-gnomAD mapping

Census counts in `data/census_canada_counts.csv` come from the 2021 Statistics Canada table 98-10-0351-01 (Canada, total age, total gender, 2021 counts); the 13 categories sum to 36,328,475. Mapping decisions are in `data/census_gnomad_mapping.csv` and implemented in `scripts/07_census_mapping.py`. Census population-group categories and gnomAD genetic ancestry groups are not equivalent, so every mapping is an approximate proxy. The mapping file was generated on 2026-10-13.

The primary mapping assigns Arab and West Asian to Middle Eastern; Black to African/African American; Chinese, Japanese and Korean to East Asian; Latin American to Admixed American; South Asian to South Asian; and "Not a visible minority" to European (non-Finnish). It covers 34,476,095 people (94.90% of the census total). Excluded are Filipino, Southeast Asian, Multiple visible minorities and Visible minority n.i.e. (1,852,380 people), because no gnomAD group is a defensible match. "Not a visible minority" includes Indigenous people and people of other ancestries, so its mapping to European (non-Finnish) is the weakest in the analysis; the selected table cannot separate Indigenous people.

Primary weights: European (non-Finnish) 0.7741, South Asian 0.0746, East Asian 0.0590, African/African American 0.0449, Middle Eastern 0.0306 and Admixed American 0.0168. They sum to 1 across six gnomAD groups. Finnish, Remaining, Amish and Ashkenazi Jewish receive no weight; their group-level results remain available separately.

Six scenarios (S0–S5) are implemented. S1 excludes "Not a visible minority" (7,786,820 people, coverage 21.43%) and is not a national estimate. S2 maps Filipino and Southeast Asian to East Asian; S3 excludes flagged groups per estimate; S4 uses exome-only frequencies; S5 uses genome-only frequencies. Finnish remains unmapped.

These estimates are descriptive approximations. Census identity categories should not be interpreted as direct measures of genetic ancestry.

## Still open

- The estimates remain exploratory population-level approximations, not clinical prescribing recommendations or estimates of the share of people who require a medication change.
- Census ancestry categories do not map directly to gnomAD ancestry groups. Indigenous people cannot be separated from the “Not a visible minority” category in the selected census table.
- European (Finnish), Remaining, Amish, and Ashkenazi Jewish gnomAD groups have no census weights in the current mapping.
- Scenario ranges are descriptive sensitivity checks, not confidence intervals or a complete quantification of mapping uncertainty.
