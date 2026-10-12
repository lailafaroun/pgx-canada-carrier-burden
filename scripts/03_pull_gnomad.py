import time
import requests
import pandas as pd
from pathlib import Path

API_URL = "https://gnomad.broadinstitute.org/api"
DATASET = "gnomad_r4"

VARIANTS = [
    {"gene": "ABCG2", "rsID": "rs2231142", "gnomad_id": "4-88131171-G-T"},
    {"gene": "VKORC1", "rsID": "rs9923231", "gnomad_id": "16-31096368-C-T"},
    {"gene": "CYP3A5", "rsID": "rs776746", "gnomad_id": "7-99672916-T-C"},
    {"gene": "CYP4F2", "rsID": "rs2108622", "gnomad_id": "19-15879621-C-T"},
    {"gene": "NUDT15", "rsID": "rs116855232", "gnomad_id": "13-48045719-C-T"},
]

POPULATIONS = {
    "afr": "African/African American",
    "amr": "Admixed American",
    "ami": "Amish",
    "asj": "Ashkenazi Jewish",
    "eas": "East Asian",
    "fin": "European (Finnish)",
    "nfe": "European (non-Finnish)",
    "mid": "Middle Eastern",
    "remaining": "Remaining",
    "sas": "South Asian",
}


def fetch_variant(variant_id):
    query = """
    query ($variantId: String!, $dataset: DatasetId!) {
      variant(variantId: $variantId, dataset: $dataset) {
        variant_id
        exome {
          ac
          an
          populations {
            id
            ac
            an
          }
        }
        genome {
          ac
          an
          populations {
            id
            ac
            an
          }
        }
      }
    }
    """

    response = requests.post(
        API_URL,
        json={
            "query": query,
            "variables": {
                "variantId": variant_id,
                "dataset": DATASET,
            },
        },
        timeout=60,
    )
    response.raise_for_status()

    result = response.json()

    if result.get("errors"):
        raise RuntimeError(f"GraphQL error for {variant_id}: {result['errors']}")

    variant = result.get("data", {}).get("variant")

    if variant is None:
        raise ValueError(f"Variant not found: {variant_id}")

    return variant


rows = []

for i, variant in enumerate(VARIANTS):
    print(f"Fetching {variant['gene']} ({variant['gnomad_id']})...")

    result = fetch_variant(variant["gnomad_id"])

    for dataset_name in ["exome", "genome"]:
        dataset_data = result.get(dataset_name)

        if not dataset_data:
            continue

        for population in dataset_data.get("populations") or []:
            pop_id = population.get("id")

            # Keep only the 10 target ancestry groups.
            # This excludes sex-specific and 1kg/HGDP subpopulation rows.
            if pop_id not in POPULATIONS:
                continue

            ac = population.get("ac")
            an = population.get("an")

            if ac is None or an is None:
                continue

            af = ac / an if an > 0 else None

            rows.append({
                "gene": variant["gene"],
                "gnomad_id": variant["gnomad_id"],
                "rsID": variant["rsID"],
                "dataset": dataset_name,
                "group": POPULATIONS[pop_id],
                "ac": ac,
                "an": an,
                "af": af,
            })

    # Pause between variants to avoid making requests too quickly.
    if i < len(VARIANTS) - 1:
        time.sleep(7)

df = pd.DataFrame(rows)

if df.empty:
    raise RuntimeError("No population data retrieved. Check the API response.")

output_path = Path("data/gnomad_variant_frequencies.csv")
output_path.parent.mkdir(parents=True, exist_ok=True)

df.to_csv(output_path, index=False)

print(f"\nSaved {len(df)} rows to {output_path}")
print("\nRows per variant and dataset:")
print(df.groupby(["gene", "dataset"]).size().to_string())

print("\nPreview:")
print(df.head(10).to_string(index=False))
