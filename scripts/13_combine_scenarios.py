from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

baseline = pd.read_csv(DATA / "national_estimates.csv")
s2 = pd.read_csv(DATA / "s2_estimates.csv")
s3 = pd.read_csv(DATA / "s3_estimates.csv")
source = pd.read_csv(DATA / "source_sensitivity_estimates.csv")

baseline["scenario"] = baseline["scenario"].map({
    "primary": "S0",
    "sensitivity": "S1",
})
assert baseline["scenario"].notna().all(), "Unexpected baseline scenario."

columns = [
    "scenario", "quantity", "national", "low", "high", "coverage",
    "weight_on_flagged_rows", "groups_excluded",
    "original_weight_excluded", "status",
]

for frame in (baseline, s2, s3, source):
    for col in columns:
        if col not in frame.columns:
            frame[col] = pd.NA

combined = pd.concat(
    [baseline[columns], s2[columns], s3[columns], source[columns]],
    ignore_index=True,
)

assert len(combined) == 48, f"Expected 48 rows, found {len(combined)}"
assert set(combined["scenario"].dropna()) == {"S0", "S1", "S2", "S3", "S4", "S5"}
assert not combined.duplicated(["scenario", "quantity"]).any()

combined = combined.sort_values(
    ["scenario", "quantity"]
).reset_index(drop=True)

combined.to_csv(DATA / "all_scenario_estimates.csv", index=False)
print(combined.groupby("scenario").size().to_string())
print(f"\nSaved {DATA / 'all_scenario_estimates.csv'}")
print("Combined scenario checks passed: 48 rows, 8 per scenario.")
