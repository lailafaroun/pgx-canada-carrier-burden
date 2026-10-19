# Project Summary

## Research question

How common are selected actionable pharmacogenomic variant alleles likely to be in
Canada, and how does the population labeling in the supporting pharmacogenomic
literature compare with selected Canadian population categories?

## Approach

The carrier analysis uses gnomAD v4.1.2 GRCh38 allele frequencies for five genes:
ABCG2, CYP3A5, CYP4F2, NUDT15, and VKORC1. Variant carrier probabilities are estimated
from allele frequencies and combined into panel estimates under an independence
assumption. Estimates are weighted using 2021 Statistics Canada population-group
counts and approximate mappings between census categories and gnomAD ancestry groups.
Sensitivity scenarios test alternative mappings, flagged low-count groups, and
exome-only or genome-only data.

A separate ClinPGx analysis counts unique publication (PMID) × population-label pairs
across variant-drug annotations for the five genes. A restricted comparison maps
selected literature labels to approximate census categories.

## Main findings

Under the primary census mapping, estimated carrier probabilities are approximately
21.7% for ABCG2, 96.3% for CYP3A5 *3, 50.1% for CYP4F2, 2.7% for NUDT15, and 60.0%
for the VKORC1 A allele. The estimated CYP3A5 expresser proportion is 21.5%, under
the model described in the project documentation. The estimated probability of
carrying at least one variant in the five-gene panel is 98.4%; this high value is
driven in part by common variants.

The evidence audit found 880 unique publication-population-label pairs. The
restricted Canada comparison includes 363 pairs across five mapped categories.
Within this subset, East Asian labels account for 68.9% of mapped evidence pairs,
compared with 26.1% of the mapped census population; Central/South Asian labels
account for 5.8% of mapped evidence pairs compared with 33.0% of the mapped census
population. These differences are descriptive and do not by themselves establish
bias or differences in evidence quality.

## Limitations and next steps

Census categories are not equivalent to genetic ancestry groups or literature
population labels. The primary mapping of “Not a visible minority” to European
(non-Finnish) is a weak proxy, and Indigenous people cannot be separately identified
from that category in the selected census table. The restricted evidence comparison
covers 21.43% of Canada's census population and is not nationally representative.
Publication-label counts are not participant counts. CYP3A5 expresser estimates
consider *3 only, and panel estimates assume independence among variants. Confidence
intervals do not capture mapping or model uncertainty.

Next steps are to review the final methods and caveats, rerun the analysis scripts,
verify generated tables and figures, and complete a final repository and Git audit.
