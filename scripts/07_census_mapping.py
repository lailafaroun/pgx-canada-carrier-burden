from pathlib import Path

import pandas as pd

COUNTS_FILE = Path("data/census_canada_counts.csv")
CARRIER_FILE = Path("data/carrier_by_group.csv")
OUTPUT_FILE = Path("data/census_gnomad_mapping.csv")

EXPECTED_TOTAL = 36_328_475
EXPECTED_PRIMARY = 34_476_095
EXPECTED_SENSITIVITY = 7_786_820

# gnomAD group, match quality, reason, included in primary analysis
DECISIONS = {
    "Arab": (
        "Middle Eastern", "approximate",
        "Census identity category used as a broad ancestry proxy.", True
    ),
    "Black": (
        "African/African American", "approximate",
        "Census identity category does not precisely represent genetic ancestry.", True
    ),
    "Chinese": (
        "East Asian", "approximate",
        "Census identity category used as a broad ancestry proxy.", True
    ),
    "Filipino": (
        "none", "unmapped",
        "No sufficiently appropriate gnomAD ancestry mapping assigned.", False
    ),
    "Japanese": (
        "East Asian", "approximate",
        "Census identity category used as a broad ancestry proxy.", True
    ),
    "Korean": (
        "East Asian", "approximate",
        "Census identity category used as a broad ancestry proxy.", True
    ),
    "Latin American": (
        "Admixed American", "approximate",
        "Broad proxy; the categories are not equivalent.", True
    ),
    "Multiple visible minorities": (
        "none", "unmapped",
        "Multiple identities cannot be assigned to a single gnomAD group.", False
    ),
    "Not a visible minority": (
        "European (non-Finnish)", "approximate",
        "Broad proxy only: this census category includes Indigenous people and others and is not equivalent to European ancestry.",
        True
    ),
    "South Asian": (
        "South Asian", "approximate",
        "Census identity category used as a broad ancestry proxy.", True
    ),
    "Southeast Asian": (
        "none", "unmapped",
        "No sufficiently appropriate gnomAD ancestry mapping assigned.", False
    ),
    "Visible minority, n.i.e.": (
        "none", "unmapped",
        "Category is not specific enough for an ancestry mapping.", False
    ),
    "West Asian": (
        "Middle Eastern", "approximate",
        "Broad proxy; West Asian and Middle Eastern are not equivalent categories.", True
    ),
}


def main():
    counts = pd.read_csv(COUNTS_FILE)
    required = {"census_group", "count"}
    assert required.issubset(counts.columns), (
        f"{COUNTS_FILE} must contain columns: {sorted(required)}"
    )
    assert counts["census_group"].is_unique, "Duplicate census groups found."
    assert set(counts["census_group"]) == set(DECISIONS), (
        "Census groups do not match the mapping decisions.\n"
        f"Missing from input: {sorted(set(DECISIONS) - set(counts['census_group']))}\n"
        f"Unexpected in input: {sorted(set(counts['census_group']) - set(DECISIONS))}"
    )

    counts["count"] = pd.to_numeric(counts["count"], errors="raise").astype(int)
    census_total = int(counts["count"].sum())
    assert census_total == EXPECTED_TOTAL, (
        f"Unexpected census total: {census_total:,}"
    )

    mapping = counts.copy()
    mapping["gnomad_group"] = mapping["census_group"].map(
        lambda x: DECISIONS[x][0]
    )
    mapping["match_quality"] = mapping["census_group"].map(
        lambda x: DECISIONS[x][1]
    )
    mapping["reason"] = mapping["census_group"].map(
        lambda x: DECISIONS[x][2]
    )
    mapping["include_primary"] = mapping["census_group"].map(
        lambda x: DECISIONS[x][3]
    )
    # Sensitivity analysis excludes "Not a visible minority".
    mapping["include_sensitivity"] = (
        mapping["include_primary"]
        & mapping["census_group"].ne("Not a visible minority")
    )

    primary = mapping[mapping["include_primary"]].copy()
    sensitivity = mapping[mapping["include_sensitivity"]].copy()

    primary_total = int(primary["count"].sum())
    sensitivity_total = int(sensitivity["count"].sum())

    assert primary_total == EXPECTED_PRIMARY, (
        f"Unexpected primary included count: {primary_total:,}"
    )
    assert sensitivity_total == EXPECTED_SENSITIVITY, (
        f"Unexpected sensitivity included count: {sensitivity_total:,}"
    )

    primary_groups = (
        primary.groupby("gnomad_group", as_index=False)["count"]
        .sum()
        .rename(columns={"count": "primary_gnomad_group_count"})
    )
    primary_groups["primary_gnomad_weight"] = (
        primary_groups["primary_gnomad_group_count"] / primary_total
    )

    sensitivity_groups = (
        sensitivity.groupby("gnomad_group", as_index=False)["count"]
        .sum()
        .rename(columns={"count": "sensitivity_gnomad_group_count"})
    )
    sensitivity_groups["sensitivity_gnomad_weight"] = (
        sensitivity_groups["sensitivity_gnomad_group_count"]
        / sensitivity_total
    )

    mapping = mapping.merge(primary_groups, on="gnomad_group", how="left")
    mapping = mapping.merge(sensitivity_groups, on="gnomad_group", how="left")

    mapping["primary_weight_contribution"] = 0.0
    mapping.loc[mapping["include_primary"], "primary_weight_contribution"] = (
        mapping.loc[mapping["include_primary"], "count"] / primary_total
    )

    mapping["sensitivity_weight_contribution"] = 0.0
    mapping.loc[mapping["include_sensitivity"], "sensitivity_weight_contribution"] = (
        mapping.loc[mapping["include_sensitivity"], "count"] / sensitivity_total
    )

    mapping["census_total"] = census_total
    mapping["primary_population_total"] = primary_total
    mapping["sensitivity_population_total"] = sensitivity_total
    mapping["primary_coverage"] = primary_total / census_total
    mapping["sensitivity_coverage"] = sensitivity_total / census_total

    # Validate that the mapped gnomAD groups exist in the carrier results.
    carrier = pd.read_csv(CARRIER_FILE)
    carrier_columns = {"group", "gnomad_group"}
    if "group" in carrier.columns:
        available_groups = set(carrier["group"].dropna().astype(str))
    elif "gnomad_group" in carrier.columns:
        available_groups = set(carrier["gnomad_group"].dropna().astype(str))
    else:
        available_groups = set()

    mapped_groups = set(primary["gnomad_group"])
    if available_groups:
        assert mapped_groups.issubset(available_groups), (
            "Mapped groups missing from carrier estimates: "
            f"{sorted(mapped_groups - available_groups)}"
        )

    assert abs(primary_groups["primary_gnomad_weight"].sum() - 1) < 1e-9
    assert abs(sensitivity_groups["sensitivity_gnomad_weight"].sum() - 1) < 1e-9
    assert mapping["reason"].notna().all()

    mapping.to_csv(OUTPUT_FILE, index=False)

    print("\nPRIMARY ANALYSIS")
    print(f"Census total: {census_total:,}")
    print(f"Included population: {primary_total:,}")
    print(f"Excluded population: {census_total - primary_total:,}")
    print(f"Mapping coverage: {primary_total / census_total:.2%}")
    print(primary_groups.to_string(index=False, formatters={
        "primary_gnomad_weight": "{:.6f}".format
    }))

    print("\nSENSITIVITY ANALYSIS (excluding 'Not a visible minority')")
    print(f"Included population: {sensitivity_total:,}")
    print(f"Excluded population: {census_total - sensitivity_total:,}")
    print(f"Mapping coverage: {sensitivity_total / census_total:.2%}")
    print(sensitivity_groups.to_string(index=False, formatters={
        "sensitivity_gnomad_weight": "{:.6f}".format
    }))

    print(f"\nSaved mapping table to {OUTPUT_FILE}")
    print("VALIDATION PASSED")


if __name__ == "__main__":
    main()
