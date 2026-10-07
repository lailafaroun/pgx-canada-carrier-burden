# Population-Level Burden of Actionable Pharmacogenomic Variants in Canada

## Overview
Clinical dosing guidelines recommend changing a drug or dose for patients who carry certain genetic variants. Variant frequencies differ substantially between ancestry groups, and most pharmacogenomic evidence comes from a limited set of populations. This project estimates how many people in Canada are likely to carry an actionable pharmacogenomic variant, and how diverse the populations are in the evidence behind the guidelines. The aim is to show where implementation and research gaps lie.

## Goal
Answer two linked questions:
1. How many people in Canada are likely to carry an actionable variant in each of 5 genes with the highest level of clinical evidence for changing prescribing (CPIC Level A), estimated from global population allele frequencies (gnomAD) and Canada's population-group composition (Statistics Canada)?
2. How diverse are the study populations in the evidence supporting these gene-drug guidelines, compared with the population the guidelines serve?

## Data Sources
- **Variant and guideline curation** - CPIC (allele tables, guidelines) and PharmGKB (gene-drug annotations, evidence)
- **Clinical significance cross-check** - ClinVar
- **Variant annotation** - Ensembl VEP
- **Allele frequencies by population group** - gnomAD, version TBD
- **Population-group counts for Canada** - Statistics Canada census table, table ID and year TBD

## Gene and Variant Selection
**Genes:** 5 genes, each with a CPIC Level A gene-drug pair (genes TBD)

**Variants:** well-defined single-nucleotide variants with rsIDs, each verified against the CPIC allele tables before inclusion.

Selected genes and variants: TBD

## Methods
1. **Curate** variants from CPIC and PharmGKB; cross-check in ClinVar
2. **Annotate** variants with Ensembl VEP
3. **Extract** allele counts and frequencies per gnomAD population group
4. **Convert** allele frequencies to carrier probabilities (Hardy-Weinberg equilibrium)
5. **Combine** variants within a gene into a gene-level carrier probability
6. **Weight** by Statistics Canada population-group counts for a national estimate, with a sensitivity analysis on group mapping
7. **Count** study populations in the evidence behind each gene-drug guideline

Data processing is in Python; figures and statistics are in R.

## Environment Setup
```bash
conda create -n pgx-canada python=3.11 pandas requests jupyter -y
conda activate pgx-canada
```
```r
install.packages(c("tidyverse"))
```

## Project Status
In progress - setup and data-availability check underway

## Repository Structure
- `README.md` - project documentation
- `scripts/` - Python data processing and R analysis scripts
- `figures/` - saved plots and visualizations
- `data/` - variant table, data-availability note, and downloaded source files

