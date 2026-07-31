"""
preprocess.py -- Cirro dataset preprocessing hook for the combined MuTect1
pipeline. Runs on the Cirro platform BEFORE Nextflow starts (via
`PreprocessDataset.from_running()`), inspects the dataset's file manifest,
and injects `mutect1_runs` into params -- the same mechanism and the same
tumor/normal discovery logic as the working Mutect2 `mutect_runs`
preprocessing script this was adapted from. Register it with this pipeline
the same way that script is registered with the Mutect2 one.

Shared-normal model (confirmed): every tumor BAM found in the dataset is
paired with the SAME single normal BAM -- there is exactly one germline/
normal sample per dataset, not one per tumor. This mirrors the Mutect2
script's active (non-commented) behavior, not its single-sample variant.

Classification heuristic: any sample whose name contains "PBMC" is treated
as the normal; everything else is treated as a tumor. This is carried over
unchanged from the Mutect2 script -- if this MuTect1 cohort's normal
samples are labeled differently, update NORMAL_NAME_MARKERS below.
"""
import json
from cirro.helpers.preprocess_dataset import PreprocessDataset

NORMAL_NAME_MARKERS = ["PBMC"]  # adjust if this cohort's normals aren't PBMC-labeled


def extract_bams(ds):
    df = ds.files.copy()
    df["file"] = df["file"].astype(str)

    df = df[df["file"].str.endswith(".bam") | df["file"].str.endswith(".bam.bai")]

    bam_map = {}

    for sample, group in df.groupby("sample"):
        bam = ""
        bai = ""

        for f in group["file"]:
            if f.endswith(".bam") and not f.endswith(".bam.bai"):
                bam = f
            elif f.endswith(".bam.bai"):
                bai = f

        if bam:
            bam_map[str(sample)] = {"bam": bam, "bai": bai}

    if not bam_map:
        raise ValueError("No BAMs found in dataset")

    return bam_map


def is_normal(sample_name):
    upper = sample_name.upper()
    return any(marker in upper for marker in NORMAL_NAME_MARKERS)


def main():
    ds = PreprocessDataset.from_running()

    print("=== ds.files preview ===")
    print(ds.files.head(20).to_string(index=False))

    bam_map = extract_bams(ds)

    normal_bam = None
    normal_bai = None
    normal_sample_name = None
    tumor_runs = []

    for sample, files in bam_map.items():
        sample_name = str(sample)

        if is_normal(sample_name):
            if normal_bam is not None:
                raise ValueError(
                    f"Multiple normal-labeled samples found ({normal_sample_name!r} and "
                    f"{sample_name!r}) -- this pipeline expects exactly one shared normal "
                    f"per dataset. Adjust NORMAL_NAME_MARKERS or the dataset contents."
                )
            normal_bam = files["bam"]
            normal_bai = files["bai"]
            normal_sample_name = sample_name
        else:
            tumor_runs.append({
                "output_prefix": sample_name,
                "tumor_reads": files["bam"],
                "tumor_reads_index": files["bai"],
                "tumor_sample_name": sample_name,
            })

    if not tumor_runs:
        raise ValueError("No tumor BAM found")

    if normal_bam is None:
        raise ValueError(
            "No normal-labeled sample found -- this pipeline requires a matched normal "
            "shared across all tumor samples (see NORMAL_NAME_MARKERS)"
        )

    mutect1_runs = []
    for run in tumor_runs:
        mutect1_runs.append({
            "output_prefix": run["output_prefix"],
            "tumor_reads": run["tumor_reads"],
            "tumor_reads_index": run["tumor_reads_index"],
            "normal_reads": normal_bam,
            "normal_reads_index": normal_bai,
            "tumor_sample_name": run["tumor_sample_name"],
            "normal_sample_name": normal_sample_name,
        })

    ds.add_param("mutect1_runs", mutect1_runs)

    print("\nFinal parameters:")
    print(json.dumps(ds.params, indent=2, default=str))


if __name__ == "__main__":
    main()