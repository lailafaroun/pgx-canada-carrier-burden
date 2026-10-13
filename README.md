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

The primary mapping assigns Canadian census groups to gnomAD ancestry groups, using European (non-Finnish) as an approximate proxy for “Not a visible minority.” The primary mapping covers 34,476,095 of 36,328,475 people (94.90%). The sensitivity scenario excludes “Not a visible minority” and covers 7,786,820 people (21.43%); it is not a national estimate.

Filipino, Southeast Asian, Multiple visible minorities, and Visible minority n.i.e. are excluded from both current scenarios. The selected census table does not separately identify Indigenous people, who are included within “Not a visible minority” along with people of other ancestries. European (Finnish), Remaining, Amish, and Ashkenazi Jewish gnomAD groups receive no census weight. Mapping Filipino and Southeast Asian groups to East Asian is planned as a Step 8 sensitivity analysis and has not yet been implemented.

These estimates are descriptive approximations. Census identity categories should not be interpreted as direct measures of genetic ancestry.

## Methods
1. **Curate** variants from CPIC and ClinPGx, cross-checked in ClinVar
2. **Annotate** variants with Ensembl VEP
3. **Extract** allele counts and frequencies per gnomAD ancestry group
4. **Convert** allele frequencies to carrier probabilities (Hardy-Weinberg equilibrium) with Wilson confidence intervals
5. **Combine** the genes into a five-gene panel probability
6. **Weight** by Statistics Canada population-group counts for a national estimate, with a sensitivity analysis on group mapping
7. **Count** study populations in the evidence behind each gene-drug guideline

Data processing is in Python; figures and statistics are in R.

## Environment Setup
```bash
conda env create -f environment.yml
conda activate pgx-canada
```
```r
renv::restore()
```

## Project Status
In progress - setup, data check, variant table, gnomAD frequency extraction, VEP annotation, carrier estimates, panel estimates and census mapping complete; weighted national estimate and sensitivity analysis next

## Repository Structure
- `README.md` - project documentation
- `scripts/` - Python data processing and R analysis scripts
- `figures/` - saved plots and visualizations
- `data/` - variant table, data-availability note, gnomAD frequency tables and downloaded source files (raw files not tracked in git)
