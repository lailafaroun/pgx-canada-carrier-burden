from pathlib import Path
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
sys.path.insert(0, str(ROOT / "scripts"))

from pgx_utils import weighted_average

WEIGHTS_FILE = DATA / "census_scenario_weights.csv"
CARRIER_FILE = DATA / "carrier_by_group.csv"
PANEL_FILE = DATA / "panel_by_group.csv"
OUTPUT_FILE = DATA / "s2_estimates.csv"

weight_df = pd.read_csv(WEIGHTS_FILE)
weight_df = weight_df[weight_df["scenario"] == "S2"].copy()
weights = weight_df[["gnomad_group", "weight"]].rename(
    columns={"gnomad_group": "group"}
)
coverage = float(weight_df["coverage"].iloc[0])

carriers = pd.read_csv(CARRIER_FILE)
panels = pd.read_csv(PANEL_FILE)

results = []

for gene in sorted(carriers["gene"].unique()):
    rows = carriers[carriers["gene"] == gene].merge(
        weights, on="group", how="inner", validate="one_to_one"
    )
    assert len(rows) == len(weights), f"Missing groups for {gene}"

    for quantity, estimate_col, low_col, high_col in [
        (gene, "carrier", "carrier_low", "carrier_high"),
    ]:
        values = rows.dropna(subset=[estimate_col, low_col, high_col])
        assert len(values) == len(weights)

        results.append({
            "quantity": quantity,
            "scenario": "S2",
            "national": weighted_average(
                values[estimate_col].tolist(), values["weight"].tolist()
            ),
            "low": weighted_average(
                values[low_col].tolist(), values["weight"].tolist()
            ),
            "high": weighted_average(
                values[high_col].tolist(), values["weight"].tolist()
            ),
            "coverage": coverage,
            "weight_on_flagged_rows": (
                values.loc[
                    values["low_an"].astype(bool) | values["few_alt"].astype(bool),
                    "weight"
                ].sum()
            ),
        })

    if gene == "CYP3A5":
        values = rows.dropna(
            subset=["expresser", "expresser_low", "expresser_high"]
        )
        results.append({
            "quantity": "CYP3A5_expressers",
            "scenario": "S2",
            "national": weighted_average(
                values["expresser"].tolist(), values["weight"].tolist()
            ),
            "low": weighted_average(
                values["expresser_low"].tolist(), values["weight"].tolist()
            ),
            "high": weighted_average(
                values["expresser_high"].tolist(), values["weight"].tolist()
            ),
            "coverage": coverage,
            "weight_on_flagged_rows": (
                values.loc[
                    values["low_an"].astype(bool) | values["few_alt"].astype(bool),
                    "weight"
                ].sum()
            ),
        })

for panel_name, label in [
    ("all_five", "all_five_panel"),
    ("without_CYP3A5_VKORC1", "reduced_panel"),
]:
    rows = panels[panels["panel"] == panel_name].merge(
        weights, on="group", how="inner", validate="one_to_one"
    )
    assert len(rows) == len(weights), f"Missing groups for {panel_name}"

    results.append({
        "quantity": label,
        "scenario": "S2",
        "national": weighted_average(
            rows["panel_carrier"].tolist(), rows["weight"].tolist()
        ),
        "low": weighted_average(
            rows["panel_low"].tolist(), rows["weight"].tolist()
        ),
        "high": weighted_average(
            rows["panel_high"].tolist(), rows["weight"].tolist()
        ),
        "coverage": coverage,
        "weight_on_flagged_rows": rows.loc[
            rows["any_low_an_or_few_alt"].astype(bool), "weight"
        ].sum(),
    })

out = pd.DataFrame(results)
assert len(out) == 8
assert (out["low"] <= out["national"]).all()
assert (out["national"] <= out["high"]).all()

out.to_csv(OUTPUT_FILE, index=False)
print(out.to_string(index=False, float_format=lambda x: f"{x:.6f}"))
print(f"\nSaved {OUTPUT_FILE.relative_to(ROOT)}")
print("S2 checks passed.")
