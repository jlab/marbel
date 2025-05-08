import shutil
import subprocess
from pathlib import Path
import os


def is_bedtools_available():
    return shutil.which("bedtools") is not None


def is_cat_available():
    return shutil.which("cat") is not None


def concat_bed_files_with_cat(dir, output_name):
    bed_files = sorted(Path(dir).glob("*.bed"))
    cmd = ["cat"] + [str(f) for f in bed_files]
    with open(output_name, "w") as out_f:
        subprocess.run(cmd, stdout=out_f, check=True)


def concat_bed_files(dir, output_name):
    bed_files = sorted(Path(dir).glob("*.bed"))
    with open(output_name, "wb") as out_f:
        for f in bed_files:
            with open(f, "rb") as in_f:
                shutil.copyfileobj(in_f, out_f)


def get_summary_paths(outdir):
    summary_dir = os.path.join(outdir, "summary")
    os.makedirs(summary_dir, exist_ok=True)
    return {
        "summary_dir": summary_dir,
        "concatted_bed": f"{summary_dir}/aggregated.bed",
        "bed": f"{summary_dir}/blocks.bed",
        "gtf": f"{summary_dir}/blocks.gtf",
        "blocks_fasta": f"{summary_dir}/blocks.fasta",
        "overlap_fasta": f"{summary_dir}/blocks_overlap.fasta",
        "overlap_tsv": f"{summary_dir}/overlap_blocks.tsv",
        "cds_ref_fasta": f"{summary_dir}/metatranscriptome_reference.fasta",
        "ref_gtf": f"{summary_dir}/metatranscriptome_reference.gtf"
    }
