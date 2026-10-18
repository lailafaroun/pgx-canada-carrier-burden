from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

all_results = pd.read_csv(DATA / "all_scenario_estimates.csv")
assert len(all_results) == 48, "Expected 48 scenario rows."

results = all_results.dropna(subset=["national"]).copy()
baseline = results[results["scenario"] == "S0"].set_index("quantity")
summary = []

for quantity, rows in results.groupby("quantity", sort=True):
    base = float(baseline.loc[quantity, "national"])

    # All available scenarios
    min_idx = rows["national"].idxmin()
    max_idx = rows["national"].idxmax()
    min_row = rows.loc[min_idx]
    max_row = rows.loc[max_idx]

    # Exclude S1, which omits the census category "Not a visible minority"
    excluding_s1 = rows[rows["scenario"] != "S1"]
    min_ex_idx = excluding_s1["national"].idxmin()
    max_ex_idx = excluding_s1["national"].idxmax()
    min_ex = excluding_s1.loc[min_ex_idx]
    max_ex = excluding_s1.loc[max_ex_idx]

    summary.append({
        "quantity": quantity,
        "S0_baseline": base,
        "minimum": float(min_row["national"]),
        "minimum_scenario": min_row["scenario"],
        "maximum": float(max_row["national"]),
        "maximum_scenario": max_row["scenario"],
        "range": float(max_row["national"] - min_row["national"]),
        "scenarios_available": int(rows["national"].count()),
        "min_excluding_S1": float(min_ex["national"]),
        "min_scenario_excluding_S1": min_ex["scenario"],
        "max_excluding_S1": float(max_ex["national"]),
        "max_scenario_excluding_S1": max_ex["scenario"],
        "range_excluding_S1": float(
            max_ex["national"] - min_ex["national"]
        ),
    })

out = pd.DataFrame(summary)
out.to_csv(DATA / "sensitivity_summary.csv", index=False)

display_cols = [
    "quantity", "S0_baseline", "minimum_scenario", "maximum_scenario",
    "range", "min_scenario_excluding_S1",
    "max_scenario_excluding_S1", "range_excluding_S1",
]
print(out[display_cols].to_string(
    index=False, float_format=lambda x: f"{x:.6f}"
))
print(f"\nSaved {DATA / 'sensitivity_summary.csv'}")
print("Sensitivity summary checks passed.")
