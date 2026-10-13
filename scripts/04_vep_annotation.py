import csv
from pathlib import Path

INPUT = Path("data/vep_results.txt")
OUTPUT = Path("data/vep_annotations.csv")
NOT_APPLICABLE = "not applicable"

VARIANTS = {
    "ABCG2": ("rs2231142", "4-88131171-G-T"),
    "VKORC1": ("rs9923231", "16-31096368-C-T"),
    "CYP3A5": ("rs776746", "7-99672916-T-C"),
    "CYP4F2": ("rs2108622", "19-15879621-C-T"),
    "NUDT15": ("rs116855232", "13-48045719-C-T"),
}

OUTPUT_FIELDS = [
    "gene", "rsID", "gnomad_id", "location", "allele", "gene_id",
    "ensembl_transcript", "refseq_mane_transcript",
    "consequence", "hgvsc", "hgvsp", "sift", "polyphen",
]

with INPUT.open(newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    rows = list(reader)

selected = []

for gene, (rsid, gnomad_id) in VARIANTS.items():
    matches = []

    for row in rows:
        extra = row.get("Extra", "")
        if "MANE_Select" not in extra:
            continue
        if rsid not in row.get("Existing_variation", ""):
            continue
        matches.append(row)

    if len(matches) != 1:
        raise ValueError(
            f"{gene} ({rsid}): expected 1 MANE Select row, "
            f"found {len(matches)}"
        )

    row = matches[0]

    # the VEP position must agree with the gnomAD ID
    chrom, pos = gnomad_id.split("-")[0], gnomad_id.split("-")[1]
    if not row.get("Location", "").startswith(f"{chrom}:{pos}"):
        raise ValueError(f"{gene}: VEP location {row.get('Location')} does not match {gnomad_id}")

    extra_fields = dict(
        item.split("=", 1)
        for item in row.get("Extra", "").split(";")
        if "=" in item
    )

    record = {
        "gene": gene,
        "rsID": rsid,
        "gnomad_id": gnomad_id,
        "location": row.get("Location", ""),
        "allele": row.get("Allele", ""),
        "gene_id": row.get("Gene", ""),
        "ensembl_transcript": row.get("Feature", ""),
        "refseq_mane_transcript": extra_fields.get("MANE_SELECT", ""),
        "consequence": row.get("Consequence", ""),
        "hgvsc": extra_fields.get("HGVSc", ""),
        "hgvsp": extra_fields.get("HGVSp", ""),
        "sift": extra_fields.get("SIFT", ""),
        "polyphen": extra_fields.get("PolyPhen", ""),
    }

    # blank means no value applies: say so instead of leaving a gap
    for key in ("hgvsc", "hgvsp", "sift", "polyphen"):
        if record[key] in ("", "-"):
            record[key] = NOT_APPLICABLE

    selected.append(record)

with OUTPUT.open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=OUTPUT_FIELDS)
    writer.writeheader()
    writer.writerows(selected)

print(f"Saved {len(selected)} annotations to {OUTPUT}")
print("Genes:", ", ".join(row["gene"] for row in selected))
