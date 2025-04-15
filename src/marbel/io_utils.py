import shutil
import subprocess
from pathlib import Path


def is_bedtools_available():
    return shutil.which("bedtools") is not None


def is_cat_available():
    return shutil.which("bedtools") is not None


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
