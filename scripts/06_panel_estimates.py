from pathlib import Path
import math

import pandas as pd


INPUT = Path("data/carrier_by_group.csv")
OUTPUT = Path("data/panel_by_group.csv")

ALL_FIVE = {"ABCG2", "CYP3A5", "CYP4F2", "NUDT15", "VKORC1"}
WITHOUT_TWO_COMMON = {"ABCG2", "CYP4F2", "NUDT15"}


def at_least_one(probabilities):
    """Combine independent probabilities of carrying at least one variant."""
    probabilities = list(probabilities)

    if any(not math.isfinite(p) or not 0 <= p <= 1
           for p in probabilities):
        raise ValueError("Probabilities must be finite and between 0 and 1.")

    return 1 - math.prod(1 - p for p in probabilities)


def run_tests():
    assert abs(at_least_one([0.1, 0.2]) - 0.28) < 1e-12
    assert at_least_one([0.37]) == 0.37
    assert at_least_one([]) == 0
    assert at_least_one([0.2, 1, 0.4]) == 1
    print("All formula tests passed.")


def flag_value(value):
    """Interpret CSV boolean values safely."""
    if pd.isna(value):
        return False
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes"}
    return bool(value)


def main():
    run_tests()

    df = pd.read_csv(INPUT)

    required = {
        "gene", "group", "carrier", "carrier_low", "carrier_high",
        "low_an", "few_alt",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    genes_found = set(df["gene"].dropna().unique())
    if genes_found != ALL_FIVE:
        raise ValueError(f"Unexpected gene set: {sorted(genes_found)}")

    numeric = ["carrier", "carrier_low", "carrier_high"]
    if df[numeric].isna().any().any():
        raise ValueError("Missing carrier estimates need investigation.")

    for col in numeric:
        if not df[col].between(0, 1).all():
            raise ValueError(f"{col} contains values outside [0, 1].")

    if not (
        (df["carrier_low"] <= df["carrier"]) &
        (df["carrier"] <= df["carrier_high"])
    ).all():
        raise ValueError("A carrier confidence interval is incorrectly ordered.")

    # One gene must contribute once per group. Duplicate records are allowed
    # only when their estimates agree; quality flags are combined with OR.
    unique_rows = []
    for (gene, group), part in df.groupby(["gene", "group"], sort=False):
        for col in numeric:
            if part[col].max() - part[col].min() > 1e-12:
                raise ValueError(
                    f"Conflicting {col} values for {gene}, {group}."
                )

        unique_rows.append({
            "gene": gene,
            "group": group,
            **{col: float(part[col].iloc[0]) for col in numeric},
            "low_an": any(part["low_an"].map(flag_value)),
            "few_alt": any(part["few_alt"].map(flag_value)),
        })

    genes = pd.DataFrame(unique_rows)

    if genes.duplicated(["gene", "group"]).any():
        raise ValueError("Gene-level input still has duplicate gene/group rows.")

    expected_groups = set(df["group"].dropna().unique())
    if len(expected_groups) != 10:
        raise ValueError(f"Expected 10 groups; found {len(expected_groups)}.")

    panel_definitions = [
        ("all_five", ALL_FIVE),
        ("without_CYP3A5_VKORC1", WITHOUT_TWO_COMMON),
    ]

    results = []

    for group in sorted(expected_groups):
        group_data = genes[genes["group"] == group]

        for panel_name, panel_genes in panel_definitions:
            included = group_data[group_data["gene"].isin(panel_genes)]

            if set(included["gene"]) != panel_genes:
                raise ValueError(
                    f"Missing gene(s) for {group}, panel {panel_name}."
                )

            results.append({
                "group": group,
                "panel": panel_name,
                "genes_included": ";".join(sorted(panel_genes)),
                "panel_carrier": at_least_one(included["carrier"]),
                "panel_low": at_least_one(included["carrier_low"]),
                "panel_high": at_least_one(included["carrier_high"]),
                "any_low_an": included["low_an"].any(),
                "any_few_alt": included["few_alt"].any(),
                "any_low_an_or_few_alt": (
                    included["low_an"] | included["few_alt"]
                ).any(),
            })

    out = pd.DataFrame(results)

    # Core output validation
    if len(out) != 20 or out.duplicated(["group", "panel"]).any():
        raise ValueError("Expected exactly 20 unique group/panel rows.")

    panel_cols = ["panel_carrier", "panel_low", "panel_high"]
    if not out[panel_cols].apply(
        lambda col: col.between(0, 1).all()
    ).all():
        raise ValueError("Panel probability outside [0, 1].")

    if not (
        (out["panel_low"] <= out["panel_carrier"]) &
        (out["panel_carrier"] <= out["panel_high"])
    ).all():
        raise ValueError("A panel confidence interval is incorrectly ordered.")

    all_five = out[out["panel"] == "all_five"].set_index("group")
    reduced = out[out["panel"] == "without_CYP3A5_VKORC1"].set_index("group")

    if not (
        reduced["panel_carrier"] <= all_five["panel_carrier"] + 1e-12
    ).all():
        raise ValueError("Reduced panel exceeds the all-five panel.")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUTPUT, index=False)

    print(f"Unique gene/group records: {len(genes)}")
    print(f"Output rows: {len(out)}")
    print(f"Saved to: {OUTPUT}")

    print("\nAll-five panel estimates:")
    print(
        all_five.reset_index()[
            ["group", "panel_carrier", "panel_low", "panel_high"]
        ].round(4).to_string(index=False)
    )

    print("\nQuality-flag counts:")
    print(out.groupby("panel")[
        ["any_low_an", "any_few_alt", "any_low_an_or_few_alt"]
    ].sum().to_string())

    print("\nAll-five panel values closest to 100%:")
    print(
        all_five.reset_index()
        .sort_values("panel_carrier", ascending=False)[
            ["group", "panel_carrier"]
        ].head(5).round(4).to_string(index=False)
    )


if __name__ == "__main__":
    main()
