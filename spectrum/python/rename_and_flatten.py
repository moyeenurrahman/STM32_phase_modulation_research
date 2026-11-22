#!/usr/bin/env python3
"""
rename_and_flatten.py
----------------------
Takes CSV files sitting inside per-case subfolders (named 1.csv, 2.csv, ...)
and copies them into the parent folder with names that analyse_spectrum.py
can match, e.g.:

    S4_/1.csv, S4_/2.csv, ...   ->   S4_1.csv, S4_2.csv, ...  (in parent folder)

Run in TWO STEPS:

    1. Leave DRY_RUN = True and run the script. It only PRINTS what it would
       do. Check the printed list carefully.

    2. Set DRY_RUN = False and run again. This actually copies the files.

Originals are left untouched - this only ever COPIES into the parent
folder, so nothing is deleted and nothing can be lost by running it wrong.
"""

import os
import re
import shutil

# ----------------------------------------------------------------------
DATA_DIR = "."          # the folder containing analyse_spectrum.py
DRY_RUN = True           # <-- set to False only after checking the preview
# ----------------------------------------------------------------------

# Folder name -> filename prefix used for matching in analyse_spectrum.py.
# A trailing underscore in the folder name is kept in the prefix, since
# analyse_spectrum.py's CASES dict expects e.g. "S3_" as the pattern.
SUBFOLDERS = [
    "S1_pointA",
    "S2_pointB",
    "N6",
    "S2N",
    "S3_",
    "S4_",
    "S5_",
    "S6N",
    "S7N",
    "S8N",
    "S8W_floor",
]


def natural_key(name):
    """Sort '2.csv' before '10.csv' instead of after."""
    return [int(t) if t.isdigit() else t
            for t in re.split(r"(\d+)", name)]


def main():
    plan = []

    for folder in SUBFOLDERS:
        folder_path = os.path.join(DATA_DIR, folder)
        if not os.path.isdir(folder_path):
            print(f"  SKIP  {folder}  (folder not found)")
            continue

        files = sorted(
            (f for f in os.listdir(folder_path) if f.lower().endswith(".csv")),
            key=natural_key,
        )
        if not files:
            print(f"  SKIP  {folder}  (no CSV files inside)")
            continue

        prefix = folder.rstrip("_")   # "S3_" -> "S3", but we re-add "_" below

        for i, fname in enumerate(files, start=1):
            src = os.path.join(folder_path, fname)
            dst_name = f"{prefix}_{i}.csv"
            dst = os.path.join(DATA_DIR, dst_name)
            plan.append((src, dst))

    if not plan:
        print("\nNothing to do - check SUBFOLDERS matches your actual folder names.")
        return

    print(f"\n{'DRY RUN - ' if DRY_RUN else ''}Planned copies ({len(plan)} files):\n")
    for src, dst in plan:
        exists = "  ALREADY EXISTS - will overwrite" if os.path.exists(dst) else ""
        print(f"  {src}  ->  {os.path.basename(dst)}{exists}")

    if DRY_RUN:
        print("\nDRY_RUN is True. Nothing was copied.")
        print("Check the list above, then set DRY_RUN = False and run again.")
        return

    print()
    for src, dst in plan:
        shutil.copy2(src, dst)
        print(f"  copied -> {os.path.basename(dst)}")

    print(f"\nDone. {len(plan)} files copied into {os.path.abspath(DATA_DIR)}")
    print("Originals in the subfolders were left untouched.")


if __name__ == "__main__":
    main()
