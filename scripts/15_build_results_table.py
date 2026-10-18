from pathlib import Path
import pandas as pd

DATA = Path("data")
OUT = Path("results/results_table.csv")
OUT.parent.mkdir(parents=True, exist_ok=True)

genes = pd.read_csv(DATA / "carrier_by_group.csv")
panels = pd.read_csv(DATA / "panel_by_group.csv")
national = pd.read_csv(DATA / "national_estimates.csv")
sensitivity = pd.read_csv(DATA / "sensitivity_summary.csv")

rows = []

# Five genes across all ten population groups.
for _, r in genes.iterrows():
    is_cyp = r["gene"] == "CYP3A5"
    rows.append({
        "gene": "CYP3A5 (expressers)" if is_cyp else r["gene"],
        "group": r["group"],
        "carrier_rate": r["expresser"] if is_cyp else r["carrier"],
        "low": r["expresser_low"] if is_cyp else r["carrier_low"],
        "high": r["expresser_high"] if is_cyp else r["carrier_high"],
        "sensitivity_low": "not applicable",
        "sensitivity_high": "not applicable",
    })

# Two panels across all ten groups.
panel_names = {
    "all_five": "all_five_panel",
    "without_CYP3A5_VKORC1": "reduced_panel",
}
for _, r in panels.iterrows():
    rows.append({
        "gene": panel_names[r["panel"]],
        "group": r["group"],
        "carrier_rate": r["panel_carrier"],
        "low": r["panel_low"],
        "high": r["panel_high"],
        "sensitivity_low": "not applicable",
        "sensitivity_high": "not applicable",
    })

# Canada weighted estimates with sensitivity ranges excluding S1.
primary = national[national["scenario"] == "primary"]
for _, r in primary.iterrows():
    s = sensitivity.loc[sensitivity["quantity"] == r["quantity"]].iloc[0]
    label = (
        "CYP3A5 (expressers)"
        if r["quantity"] == "CYP3A5_expressers"
        else r["quantity"]
    )
    rows.append({
        "gene": label,
        "group": "Canada (weighted)",
        "carrier_rate": r["national"],
        "low": r["low"],
        "high": r["high"],
        "sensitivity_low": s["min_excluding_S1"],
        "sensitivity_high": s["max_excluding_S1"],
    })

result = pd.DataFrame(rows)

# Checks: ten groups, seven quantities per group, and eight Canadian rows.
expected = {
    "ABCG2", "CYP3A5 (expressers)", "CYP4F2", "NUDT15", "VKORC1",
    "all_five_panel", "reduced_panel",
}
groups = set(genes["group"].unique())
assert len(groups) == 10
for group in groups:
    assert set(result.loc[result["group"] == group, "gene"]) == expected

canada = result[result["group"] == "Canada (weighted)"]
assert len(canada) == 8, f"Expected 8 Canadian rows; found {len(canada)}"
assert set(canada["gene"]) == expected | {"CYP3A5"}

assert not result[["carrier_rate", "low", "high"]].isna().any().any()
assert not result[["sensitivity_low", "sensitivity_high"]].isna().any().any()

result = result[
    ["gene", "group", "carrier_rate", "low", "high",
     "sensitivity_low", "sensitivity_high"]
]
result.to_csv(OUT, index=False)

print(f"Saved: {OUT}")
print(f"Total rows: {len(result)}")
print(f"Population groups: {len(groups)}")
print(f"Canada weighted rows: {len(canada)}")
print("All checks passed.")
print("\nCanada weighted estimates:")
print(canada.to_string(index=False))
