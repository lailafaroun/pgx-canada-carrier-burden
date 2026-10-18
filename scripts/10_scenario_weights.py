from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

COUNTS_FILE = DATA / "census_canada_counts.csv"
MAPPING_FILE = DATA / "census_mapping_scenarios.csv"
OUTPUT_FILE = DATA / "census_scenario_weights.csv"

counts = pd.read_csv(COUNTS_FILE)
mapping = pd.read_csv(MAPPING_FILE)

counts["count"] = pd.to_numeric(counts["count"], errors="raise")

assert counts["census_group"].is_unique, "Duplicate census groups found."
assert mapping[["scenario", "census_group"]].duplicated().sum() == 0, (
    "Duplicate scenario/census mappings found."
)

df = mapping.merge(counts, on="census_group", how="left", validate="many_to_one")
assert df["count"].notna().all(), "Some mapping categories lack census counts."

results = []

for scenario, rows in df.groupby("scenario", sort=True):
    mapped = rows.dropna(subset=["gnomad_group"]).copy()
    mapped["gnomad_group"] = mapped["gnomad_group"].astype(str)

    total_census = int(rows["count"].sum())
    mapped_total = int(mapped["count"].sum())

    assert mapped_total > 0, f"{scenario} has no mapped population."

    grouped = (
        mapped.groupby("gnomad_group", as_index=False)["count"]
        .sum()
        .rename(columns={"count": "group_count"})
    )
    grouped["scenario"] = scenario
    grouped["weight"] = grouped["group_count"] / mapped_total
    grouped["coverage"] = mapped_total / total_census

    assert abs(grouped["weight"].sum() - 1.0) < 1e-9

    results.append(grouped[
        ["scenario", "gnomad_group", "group_count", "weight", "coverage"]
    ])

out = pd.concat(results, ignore_index=True)
out.to_csv(OUTPUT_FILE, index=False)

print(out.to_string(index=False, float_format=lambda x: f"{x:.6f}"))
print(f"\nSaved {OUTPUT_FILE.relative_to(ROOT)}")
print("Scenario weight checks passed.")
