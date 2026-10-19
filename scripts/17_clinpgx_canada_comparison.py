from pathlib import Path
import pandas as pd

DATA = Path("data")
evidence = pd.read_csv(DATA / "evidence_counts.csv")
census = pd.read_csv(DATA / "census_canada_counts.csv")

# Conservative crosswalk: approximate identity/population-label comparisons,
# not claims that census identity categories equal genetic ancestry.
crosswalk = [
    ("European", "", "Excluded", "Not a visible minority is not equivalent to European ancestry."),
    ("East Asian", "Chinese;Japanese;Korean", "Approximate", "Census categories do not fully represent all East Asian populations."),
    ("African American/Afro-Caribbean", "Black", "Approximate", "Broad population-label proxy; categories are not equivalent."),
    ("Sub-Saharan African", "Black", "Approximate", "Broad population-label proxy; aggregated with the other Black-mapped label."),
    ("Near Eastern", "Arab;West Asian", "Approximate", "Broad regional proxy; categories are not equivalent."),
    ("Central/South Asian", "South Asian", "Approximate", "Broad regional proxy; categories are not equivalent."),
    ("Latino", "Latin American", "Approximate", "Broad proxy; labels and census identity categories differ."),
    ("American", "", "Excluded", "Label is ambiguous and cannot be assigned defensibly."),
    ("Multiple groups", "", "Excluded", "Cannot assign to a single census category."),
    ("Unknown", "", "Excluded", "Population label is unknown or unspecified."),
    ("Non-population / in vitro", "", "Excluded", "Not a population category."),
]

mapping = pd.DataFrame(
    crosswalk,
    columns=["clinpgx_group", "census_group", "mapping_status", "reason"]
)
mapping["include_in_comparison"] = mapping["census_group"].ne("")
mapping.to_csv(DATA / "clinpgx_canada_mapping.csv", index=False)

# Recompute the total and validate the evidence-share denominator.
total_pairs = int(evidence["unique_pmid_group_pairs"].sum())
assert total_pairs == 880, f"Unexpected total PMID-label pairs: {total_pairs}"

# Independent hand check for European label.
european = int(
    evidence.loc[evidence["group"].eq("European"),
                 "unique_pmid_group_pairs"].iloc[0]
)
hand_share = european / total_pairs
assert abs(hand_share - 126 / 880) < 1e-12
print(f"Hand check: European = {european}/{total_pairs} = {hand_share:.4%}")

# Aggregate evidence labels that map to the same census category.
census_groups = {
    "Chinese": "East Asian",
    "Japanese": "East Asian",
    "Korean": "East Asian",
    "Black": "Black",
    "Arab": "Near Eastern",
    "West Asian": "Near Eastern",
    "South Asian": "Central/South Asian",
    "Latin American": "Latino",
}
census["comparison_group"] = census["census_group"].map(census_groups)
mapped_census = census.dropna(subset=["comparison_group"]).copy()

ev_map = mapping[mapping["include_in_comparison"]].copy()
ev_map["census_group"] = ev_map["census_group"].str.split(";")
ev_map = ev_map.explode("census_group")
ev_map["comparison_group"] = ev_map["census_group"].map({
    "Chinese": "East Asian",
    "Japanese": "East Asian",
    "Korean": "East Asian",
    "Black": "Black",
    "Arab": "Near Eastern",
    "West Asian": "Near Eastern",
    "South Asian": "Central/South Asian",
    "Latin American": "Latino",
})

# One evidence label may map to multiple census categories, so aggregate
# those census counts first, and use one evidence total per comparison group.
evidence_by_group = (
    evidence.merge(
        ev_map[["clinpgx_group", "comparison_group"]].drop_duplicates(),
        left_on="group", right_on="clinpgx_group", how="inner"
    )
    .groupby("comparison_group", as_index=False)["unique_pmid_group_pairs"]
    .sum()
    .rename(columns={"unique_pmid_group_pairs": "evidence_pmid_group_pairs"})
)

census_by_group = (
    mapped_census.groupby("comparison_group", as_index=False)["count"]
    .sum()
    .rename(columns={"count": "census_population_count"})
)

comparison = evidence_by_group.merge(census_by_group, on="comparison_group")
comparison["evidence_share_mapped"] = (
    comparison["evidence_pmid_group_pairs"]
    / comparison["evidence_pmid_group_pairs"].sum()
)
comparison["census_share_mapped"] = (
    comparison["census_population_count"]
    / comparison["census_population_count"].sum()
)
comparison["evidence_minus_census_percentage_points"] = (
    comparison["evidence_share_mapped"] - comparison["census_share_mapped"]
) * 100

comparison.to_csv(DATA / "evidence_canada_comparison.csv", index=False)

print(f"\nTotal evidence PMID-label pairs: {total_pairs}")
print(f"Unknown share: {181 / total_pairs:.2%}")
print("\nCrosswalk:")
print(mapping.to_string(index=False))
print("\nEvidence vs. census (shares calculated within mapped categories only):")
print(comparison.round(4).to_string(index=False))
print("\nSaved:")
print("  data/clinpgx_canada_mapping.csv")
print("  data/evidence_canada_comparison.csv")
