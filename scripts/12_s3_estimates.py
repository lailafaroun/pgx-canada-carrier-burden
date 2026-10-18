from pathlib import Path
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
sys.path.insert(0, str(ROOT / "scripts"))

from pgx_utils import weighted_average

weights_df = pd.read_csv(DATA / "census_scenario_weights.csv")
weights_df = weights_df[weights_df["scenario"] == "S0"].copy()
weights = weights_df[["gnomad_group", "weight"]].rename(
    columns={"gnomad_group": "group"}
)
base_coverage = float(weights_df["coverage"].iloc[0])

carriers = pd.read_csv(DATA / "carrier_by_group.csv")
panels = pd.read_csv(DATA / "panel_by_group.csv")
results = []

def calculate(quantity, rows, estimate_col, low_col, high_col, flagged):
    original_weight_excluded = float(rows.loc[flagged, "weight"].sum())
    kept = rows.loc[~flagged].copy()

    assert not kept.empty, f"No unflagged groups remain for {quantity}"
    kept["weight"] = kept["weight"] / kept["weight"].sum()

    national = weighted_average(
        kept[estimate_col].tolist(), kept["weight"].tolist()
    )
    low = weighted_average(kept[low_col].tolist(), kept["weight"].tolist())
    high = weighted_average(kept[high_col].tolist(), kept["weight"].tolist())

    results.append({
        "quantity": quantity,
        "scenario": "S3",
        "national": national,
        "low": low,
        "high": high,
        "coverage": base_coverage * (1 - original_weight_excluded),
        "weight_on_flagged_rows": 0.0,
        "groups_excluded": int(flagged.sum()),
        "original_weight_excluded": original_weight_excluded,
    })

for gene in sorted(carriers["gene"].unique()):
    rows = carriers[carriers["gene"] == gene].merge(
        weights, on="group", how="inner", validate="one_to_one"
    )
    assert len(rows) == len(weights), f"Missing groups for {gene}"

    flagged = rows["low_an"].astype(bool) | rows["few_alt"].astype(bool)
    calculate(gene, rows, "carrier", "carrier_low", "carrier_high", flagged)

    if gene == "CYP3A5":
        calculate(
            "CYP3A5_expressers", rows,
            "expresser", "expresser_low", "expresser_high", flagged
        )

for panel_name, label in [
    ("all_five", "all_five_panel"),
    ("without_CYP3A5_VKORC1", "reduced_panel"),
]:
    rows = panels[panels["panel"] == panel_name].merge(
        weights, on="group", how="inner", validate="one_to_one"
    )
    assert len(rows) == len(weights), f"Missing groups for {panel_name}"

    flagged = rows["any_low_an_or_few_alt"].astype(bool)
    calculate(label, rows, "panel_carrier", "panel_low", "panel_high", flagged)

out = pd.DataFrame(results)
assert len(out) == 8
assert (out["low"] <= out["national"]).all()
assert (out["national"] <= out["high"]).all()

out.to_csv(DATA / "s3_estimates.csv", index=False)
print(out.to_string(index=False, float_format=lambda x: f"{x:.6f}"))
print(f"\nSaved {DATA / 's3_estimates.csv'}")
print("S3 checks passed.")
