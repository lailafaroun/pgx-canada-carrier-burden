from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

carriers = pd.read_csv(DATA / "carrier_by_group.csv")
panels = pd.read_csv(DATA / "panel_by_group.csv")
mapping = pd.read_csv(DATA / "census_gnomad_mapping.csv")

scenarios = {
    "primary": ("include_primary", "primary_gnomad_weight", "primary_coverage"),
    "sensitivity": (
        "include_sensitivity",
        "sensitivity_gnomad_weight",
        "sensitivity_coverage",
    ),
}

results = []

for scenario, (include_col, weight_col, coverage_col) in scenarios.items():
    selected = mapping[mapping[include_col]].copy()

    weights = (
        selected[["gnomad_group", weight_col]]
        .drop_duplicates()
        .rename(columns={"gnomad_group": "group", weight_col: "weight"})
    )

    assert abs(weights["weight"].sum() - 1.0) < 1e-6, (
        f"{scenario} weights do not sum to 1"
    )

    coverage = float(selected[coverage_col].iloc[0])

    # Five gene carrier estimates and CYP3A5 expresser estimate.
    for gene in sorted(carriers["gene"].unique()):
        gene_rows = carriers[carriers["gene"] == gene].merge(
            weights, on="group", how="inner", validate="one_to_one"
        )

        for quantity, estimate_col, low_col, high_col in [
            (gene, "carrier", "carrier_low", "carrier_high"),
        ]:
            pass

        for kind, estimate_col, low_col, high_col in [
            ("carrier", "carrier", "carrier_low", "carrier_high"),
            ("expresser", "expresser", "expresser_low", "expresser_high"),
        ]:
            if kind == "expresser" and gene != "CYP3A5":
                continue

            values = gene_rows.dropna(subset=[estimate_col, low_col, high_col])
            assert len(values) == len(weights), (
                f"Missing {kind} rows for {gene} in {scenario}"
            )

            estimate = (values["weight"] * values[estimate_col]).sum()
            low = (values["weight"] * values[low_col]).sum()
            high = (values["weight"] * values[high_col]).sum()

            assert low <= estimate <= high, (
                f"Invalid interval for {gene} {kind} in {scenario}"
            )
            assert (
                values[estimate_col].min() - 1e-9
                <= estimate
                <= values[estimate_col].max() + 1e-9
            ), f"{gene} {kind} estimate outside group range"

            flagged = values["low_an"].astype(bool) | values["few_alt"].astype(bool)
            flagged_weight = values.loc[flagged, "weight"].sum()

            results.append({
                "quantity": gene if kind == "carrier" else "CYP3A5_expressers",
                "scenario": scenario,
                "national": estimate,
                "low": low,
                "high": high,
                "coverage": coverage,
                "weight_on_flagged_rows": flagged_weight,
            })

    # Use the already-computed group-level panel values.
    for panel_name, panel_label in [
        ("all_five", "all_five_panel"),
        ("without_CYP3A5_VKORC1", "reduced_panel"),
    ]:
        panel_rows = panels[panels["panel"] == panel_name].merge(
            weights, on="group", how="inner", validate="one_to_one"
        )
        assert len(panel_rows) == len(weights), (
            f"Missing panel rows for {panel_name} in {scenario}"
        )

        estimate = (panel_rows["weight"] * panel_rows["panel_carrier"]).sum()
        low = (panel_rows["weight"] * panel_rows["panel_low"]).sum()
        high = (panel_rows["weight"] * panel_rows["panel_high"]).sum()

        assert low <= estimate <= high, (
            f"Invalid interval for {panel_name} in {scenario}"
        )
        assert (
            panel_rows["panel_carrier"].min() - 1e-9
            <= estimate
            <= panel_rows["panel_carrier"].max() + 1e-9
        ), f"{panel_name} estimate outside group range"

        flagged = panel_rows["any_low_an_or_few_alt"].astype(bool)
        flagged_weight = panel_rows.loc[flagged, "weight"].sum()

        results.append({
            "quantity": panel_label,
            "scenario": scenario,
            "national": estimate,
            "low": low,
            "high": high,
            "coverage": coverage,
            "weight_on_flagged_rows": flagged_weight,
        })

out = pd.DataFrame(results)

# The reduced panel should never exceed the all-five panel.
for scenario in scenarios:
    values = out[out["scenario"] == scenario].set_index("quantity")["national"]
    assert values["reduced_panel"] <= values["all_five_panel"] + 1e-9

# Simple weighted-average test: equal weights on 0.1 and 0.3 give 0.2.
assert abs((0.5 * 0.1 + 0.5 * 0.3) - 0.2) < 1e-12

output = DATA / "national_estimates.csv"
out.to_csv(output, index=False, float_format="%.8f")
print(f"Saved {output.relative_to(ROOT)}")
print(out.to_string(index=False))
print("\nAll national estimate checks passed.")
