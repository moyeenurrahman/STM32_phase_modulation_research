"""
stage.py
========
Builds a ready-to-run work/ folder for the phase accuracy analysis scripts.

Why this is needed
------------------
The scripts were written while the data sat in one flat folder, with
lower-case ".csv" names such as "m1_kp1p00_r1.csv" and "ref_carrier_kp0.csv".
In this repository the CSV files are sorted into sub-folders and have
upper-case ".CSV" extensions (Linux and macOS treat these as different
files). This script copies everything into one flat folder with the names the
scripts expect. Nothing in the repository is changed.

Usage (from the phase_accuracy folder)
--------------------------------------
    python stage.py
    cd work
    python measure_fc_ref.py ref_carrier_kp0.csv
    ...                       (see phase_accuracy_procedure.md, Section 6)
"""

import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC_PY = os.path.join(HERE, "python")
SRC_CSV = os.path.join(HERE, "csv")
WORK = os.path.join(HERE, "work")

KP_M1 = ["0.25", "0.50", "0.75", "1.00", "1.50", "2.00", "2.50", "3.00",
         "3.50", "4.00"]

missing = []


def copy(src, dst_name):
    if not os.path.isfile(src):
        missing.append(os.path.relpath(src, HERE))
        return False
    shutil.copyfile(src, os.path.join(WORK, dst_name))
    return True


def main():
    if os.path.isdir(WORK):
        shutil.rmtree(WORK)
    os.makedirs(WORK)

    # ---- scripts (everything in python/, manifests are rewritten below) ----
    n_py = 0
    for name in sorted(os.listdir(SRC_PY)):
        if name.endswith(".py"):
            shutil.copyfile(os.path.join(SRC_PY, name),
                            os.path.join(WORK, name))
            n_py += 1
    for need in ("measure_fc_ref.py", "check_m1_point.py",
                 "check_m3_point.py", "diagnose.py", "check_r2.py",
                 "check_r3_noise.py", "check_m2_point.py"):
        if not os.path.isfile(os.path.join(WORK, need)):
            missing.append("python/" + need)

    # ---- reference captures ----
    ref = os.path.join(SRC_CSV, "Reference")
    copy(os.path.join(ref, "R1.CSV"), "ref_carrier_kp0.csv")
    for n in ("R2_CH1.CSV", "R2_CH2.CSV", "R3.CSV", "R4.CSV"):
        copy(os.path.join(ref, n), n)

    # ---- PH1 (M1): 10 kp values x 3 repeats ----
    m1_rows = []
    for kp in KP_M1:
        tag = kp.replace(".", "p")
        for r in (1, 2, 3):
            name = "m1_kp%s_r%d.csv" % (tag, r)
            if copy(os.path.join(SRC_CSV, "M1",
                                 "m1_kp%s_r%d.CSV" % (tag, r)), name):
                m1_rows.append("%s  %s" % (kp, name))
    with open(os.path.join(WORK, "m1_manifest.txt"), "w") as fh:
        fh.write('# PH1 (M1) manifest, one capture per line: "kp  filename"\n')
        fh.write("# Written by stage.py\n")
        fh.write("\n".join(m1_rows) + "\n")

    # ---- PH2 (M2): coarse and fine ----
    for kp in (2.20, 2.25, 2.30, 2.35, 2.40, 2.45, 2.50, 2.55, 2.60):
        copy(os.path.join(SRC_CSV, "M2", "J0 coarse pass",
                          "m2_j0_coarse_%.2f.CSV" % kp),
             "m2_j0_coarse_%.2f.csv" % kp)
    for kp in (2.30, 2.32, 2.34, 2.36, 2.38, 2.40, 2.42, 2.44, 2.46, 2.48,
               2.50):
        copy(os.path.join(SRC_CSV, "M2", "J0 fine pass",
                          "m2_j0_fine_%.2f.CSV" % kp),
             "m2_j0_fine_%.2f.csv" % kp)

    # ---- PH3 (M3): three repeats and the zoomed transition ----
    m3_rows = []
    for r in (1, 2, 3):
        a, b = "m3_ch1_r%d.csv" % r, "m3_ch2_r%d.csv" % r
        ok1 = copy(os.path.join(SRC_CSV, "M3", "m3_ch1_r%d.CSV" % r), a)
        ok2 = copy(os.path.join(SRC_CSV, "M3", "m3_ch2_r%d.CSV" % r), b)
        if ok1 and ok2:
            m3_rows.append("%s  %s" % (a, b))
    copy(os.path.join(SRC_CSV, "M3", "m3_transition.CSV"),
         "m3_transition.csv")
    copy(os.path.join(SRC_CSV, "M3", "m3_transition_ch2.CSV"),
         "m3_transition_ch2.csv")
    with open(os.path.join(WORK, "m3_manifest.txt"), "w") as fh:
        fh.write('# PH3 (M3) manifest, one repeat per line: '
                 '"ch1_file  ch2_file"\n# Written by stage.py\n')
        fh.write("\n".join(m3_rows) + "\n")

    n_csv = len([f for f in os.listdir(WORK) if f.lower().endswith(".csv")])
    print("Staged %d scripts and %d CSV files into %s" % (n_py, n_csv, WORK))
    if missing:
        print("\nMISSING (not found, not staged):")
        for m in missing:
            print("  ", m)
        sys.exit(1)
    print("Done. Now:  cd work")


if __name__ == "__main__":
    main()
