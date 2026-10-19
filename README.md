# Population-Level Burden of Actionable Pharmacogenomic Variants in Canada

## Overview
Clinical dosing guidelines recommend changing a drug or dose for patients who carry certain genetic variants. Variant frequencies differ substantially between ancestry groups, and most pharmacogenomic evidence comes from a limited set of populations. This project estimates how many people in Canada are likely to carry an actionable pharmacogenomic variant, and how diverse the populations are in the evidence behind the guidelines. The aim is to show where implementation and research gaps lie.

## Goal
Answer two linked questions:
1. How many people in Canada are likely to carry an actionable variant in each of 5 genes with the highest level of clinical evidence for changing prescribing (CPIC Level A), estimated from global population allele frequencies (gnomAD) and Canada's population-group composition (Statistics Canada)?
2. How diverse are the study populations in the evidence supporting these gene-drug guidelines, compared with the population the guidelines serve?

## Data Sources
- **Variant and guideline curation** - CPIC (allele tables, guidelines) and ClinPGx, formerly PharmGKB (gene-drug annotations, evidence)
- **Clinical significance cross-check** - ClinVar
- **Variant annotation** - Ensembl VEP
- **Allele frequencies by population group** - gnomAD v4.1.2 (GRCh38)
- **Population-group counts for Canada** - Statistics Canada 2021 census, table 98-10-0351-01

## Gene and Variant Selection
**Genes:** 5 genes, each with a CPIC Level A gene-drug pair.

**Variants:** one well-defined single-nucleotide variant with an rsID per gene, verified against the CPIC allele tables, present in gnomAD v4, and common enough in at least one ancestry group to estimate.

| Gene | Variant | rsID | CPIC Level A drug(s) |
|---|---|---|---|
| ABCG2 | Q141K | rs2231142 | rosuvastatin |
| VKORC1 | -1639G>A | rs9923231 | warfarin |
| CYP3A5 | `*3` | rs776746 | tacrolimus |
| CYP4F2 | `*3` (V433M) | rs2108622 | warfarin |
| NUDT15 | `*3` (R139C) | rs116855232 | azathioprine, mercaptopurine, thioguanine |

## Analysis Decisions
1. **gnomAD data:** exomes and genomes combined (AC and AN added) as the primary dataset; exomes only (four genes) and genomes only (all five genes) as sensitivity checks. The VKORC1 variant is in genomes only.
2. **What is counted:** the share of people carrying at least one copy of the variant allele (gnomAD's alternate allele), 1 − (1 − AF)². For CYP3A5, a second estimate gives expressers (at least one working copy), 1 − AF(`*3`)². These are carrier estimates, not the share needing a prescribing change. The five-gene panel is shown with and without the two most common variants.
3. **Groups:** all 10 gnomAD genetic ancestry groups are kept with original labels; small samples are flagged, not excluded.
4. **Precision flags:** rows (variant × group) with an allele number below 2,500 are flagged `low_an` (about ±2 percentage points at a frequency of 50%); rows with fewer than 20 alternate alleles are flagged `few_alt`; AC = 0 means not observed. Flagged rows are kept and marked.
5. **Panel:** the probability of carrying at least one of the five variants is 1 − ∏(1 − cᵢ), assuming independence (the variants are on different chromosomes, but population structure can still correlate them). It is reported for all five genes and without CYP3A5 and VKORC1, the two with the highest median carrier probability; the all-five panel is at or above 99% in eight of ten groups. The panel is carriage of at least one variant, not the share needing a prescribing change.

## Census-to-gnomAD mapping

The primary mapping assigns Canadian census groups to gnomAD ancestry groups, using European (non-Finnish) as an approximate proxy for “Not a visible minority.” The primary mapping covers 34,476,095 of 36,328,475 people (94.90%). S1 excludes “Not a visible minority” and covers 7,786,820 people (21.43%); it describes only the mapped subset and is not a national estimate.

S2 assigns Filipino and Southeast Asian census groups to East Asian, increasing mapped coverage to 98.61%. S3 excludes flagged gnomAD groups separately for each estimate and renormalizes the remaining weights. S4 uses exome-only frequencies and S5 uses genome-only frequencies; VKORC1 is unavailable in exomes. These scenarios test mapping, flagged data, and sequencing-source choices.

The selected census table does not separately identify Indigenous people, who are included within “Not a visible minority” along with people of other ancestries. European (Finnish), Remaining, Amish, and Ashkenazi Jewish gnomAD groups receive no census weight.

## Methods
1. **Curate** variants from CPIC and ClinPGx, cross-checked in ClinVar
2. **Annotate** variants with Ensembl VEP
3. **Extract** allele counts and frequencies per gnomAD ancestry group
4. **Convert** allele frequencies to carrier probabilities (Hardy-Weinberg equilibrium) with Wilson confidence intervals
5. **Combine** the genes into a five-gene panel probability
6. **Weight** by Statistics Canada population-group counts for a national estimate, with a sensitivity analysis on group mapping
7. **Count** study populations in the evidence behind each gene-drug guideline

The consolidated results table combines gene-level carrier estimates, CYP3A5 expresser estimates, panel estimates, and weighted Canada estimates. Canada rows include the sensitivity range excluding S1, which omits the largest census category. Data processing is in Python; figures and statistics are in R.

## Environment Setup
```bash
conda env create -f environment.yml
conda activate pgx-canada
```
```r
renv::restore()
```

## Results: National Estimates

Weighted national estimates of the share of people carrying at least one copy of each variant allele (gnomAD alternate allele), under the primary census mapping, with additional mapping, flagged-data, and sequencing-source sensitivity scenarios reported separately. Intervals combine the group-level Wilson bounds, so they reflect gnomAD sampling uncertainty only, not census, mapping or model assumptions.

| Quantity | Primary mapping (94.9% coverage) | S1: excludes “Not a visible minority” (21.4% coverage) |
|---|---|---|
| ABCG2 (Q141K carriers) | 21.7% (21.5-21.9) | 24.1% (23.6-24.6) |
| VKORC1 (-1639G>A carriers) | 60.0% (59.3-60.7) | 55.6% (54.0-57.1) |
| CYP3A5 (`*3` carriers) | 96.3% (96.3-96.4) | 85.3% (85.0-85.6) |
| CYP3A5 expressers (at least one working copy) | 21.5% (21.3-21.8) | 50.9% (50.2-51.6) |
| CYP4F2 (V433M carriers) | 50.1% (49.8-50.3) | 48.0% (47.3-48.6) |
| NUDT15 (`*3` carriers) | 2.7% (2.6-2.8) | 10.3% (10.0-10.7) |
| At least one of the five variants | 98.4% (98.4-98.5) | 93.3% (93.1-93.5) |
| At least one of ABCG2, CYP4F2, NUDT15 | 61.7% (61.4-62.0) | 63.2% (62.5-64.0) |

- The estimates depend strongly on the weights. The primary mapping is about 77% European (non-Finnish), so variants that are rare in Europeans look rare nationally: NUDT15 is 2.7% (primary) against 10.3% under S1, which excludes “Not a visible minority”, and CYP3A5 expressers are 21.5% against 50.9%.
- The five-gene panel is 98.4% in the primary mapping, driven by the very common CYP3A5 `*3` and VKORC1 variants, so it separates little; the panel without them is 61.7%.
- VKORC1 and the all-five panel put 3.1% (primary) and 13.5% (S1: excludes “Not a visible minority”) of their weight on a group with a flagged row (VKORC1 in Middle Eastern, low allele number).
- These are carriage estimates, not the share needing a prescribing change. In the primary mapping, "Not a visible minority" is assigned to European (non-Finnish) and includes Indigenous people; S1, which excludes “Not a visible minority,” is not a national estimate.

## Pharmacogenomic evidence diversity

ClinPGx study-population labels were summarized as unique publication (PMID) ×
population-label pairs across annotations for the five genes. A publication may
contribute to more than one population label. The analysis includes records labeled
Unknown, Multiple groups, and Non-population / in vitro; these categories are not
equivalent to identified ancestry groups.

The Canada comparison maps only selected ClinPGx labels to approximate Statistics
Canada population categories. European evidence is excluded from this comparison
because it cannot be defensibly matched to the selected census categories; ambiguous
and unclassified evidence labels are also excluded. The comparison contains 363
mapped publication-population pairs out of 880 total pairs. The corresponding census
categories cover 7,786,820 people, or 21.43% of Canada's 2021 population in the
selected table. Evidence shares and census shares are each normalized within their
own mapped subset.

This is a descriptive comparison of evidence labels, not a count of study participants
or a direct measure of research bias, evidence quality, or clinical outcomes. ClinPGx
population labels and census identity categories are not equivalent. The restricted
comparison is not nationally representative. See
`data/evidence_canada_comparison.csv` and `data/clinpgx_canada_mapping.csv`.

## Results: Sensitivity Analysis

S0 is the primary mapping. S1 excludes “Not a visible minority” and has only 21.43% coverage, so its range is shown separately from the range across S0 and S2–S5. These are scenario ranges, not confidence intervals.

| Quantity | S0 baseline | Range, all scenarios | Range excluding S1 |
|---|---:|---:|---:|
| ABCG2 | 21.669% | 2.663% | 1.318% |
| CYP3A5 | 96.336% | 11.015% | 0.122% |
| CYP3A5_expressers | 21.540% | 29.443% | 0.955% |
| CYP4F2 | 50.068% | 2.175% | 1.296% |
| NUDT15 | 2.731% | 7.616% | 0.634% |
| VKORC1 | 60.009% | 5.891% | 2.038% |
| all_five_panel | 98.431% | 5.206% | 0.106% |
| reduced_panel | 61.666% | 2.560% | 1.668% |

S2 tests an expanded census mapping; S3 excludes flagged population groups per estimate; S4 and S5 compare exome-only and genome-only frequencies. The ranges describe sensitivity to these analysis choices and should not be interpreted as sampling uncertainty.
## Results Table

`results/results_table.csv` has 88 rows: 8 quantities for each of the 10 gnomAD groups, plus 8 rows for Canada (weighted, primary mapping).

| Column | Meaning |
|---|---|
| `quantity` | ABCG2, VKORC1, CYP3A5, CYP4F2 and NUDT15: carrying at least one copy of the variant allele; `CYP3A5_expressers`: at least one working copy; `all_five_panel` and `reduced_panel`: at least one of the variants |
| `group` | gnomAD genetic ancestry group, or `Canada (weighted)` |
| `estimate`, `low`, `high` | Estimate with 95% interval (Wilson bounds from the allele counts; gnomAD sampling uncertainty only) |
| `an` | Allele number for the gene rows (blank for panels and Canada) |
| `flagged` | True if the row has an allele number below 2,500 or fewer than 20 alternate alleles; for panels, if any contributing row is flagged |
| `weight_on_flagged_rows`, `coverage` | Canada rows: share of the census weight on flagged groups; share of the census population covered by the mapping |
| `sensitivity_low`, `sensitivity_high` | Canada rows: smallest and largest national estimate across the sensitivity scenarios other than S1 |

## Project Status
Analysis complete: variant curation, gnomAD carrier estimates, panel estimates, Canadian census weighting, sensitivity analysis, ClinPGx evidence counts, evidence-to-census comparison, and figures. See SUMMARY.md for the main findings and limitations.

## Repository Structure
- `README.md` - project documentation
- `scripts/` - Python data processing and R analysis scripts
- `figures/` - saved plots and visualizations
- `results/` - results table: estimates by gene and gnomAD group and for Canada (weighted), with intervals and sensitivity ranges
- `data/` - variant table, data-availability note, gnomAD frequency tables, census and panel tables, and downloaded source files (raw files not tracked in git)
