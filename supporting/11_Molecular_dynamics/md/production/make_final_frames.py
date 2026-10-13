#!/usr/bin/env python3
"""The last solute frame of each of the eleven long runs, for git and for the supplementary archive (md/PLAN.md,
section 6: "inputs, parameters, final frames and analysis tables go into git").

The driver (md/run_md.py) wrote `final_solute.pdb` into each run folder at the end of the last production
segment: the 4,094 atoms of the solute (protein with its acetyl cap, zinc, ligand) at 100 ns, coordinates in
Angstrom as integrated.  This script compresses each file with gzip, unchanged, into md/results/final_frames/.
The uncompressed content is byte for byte the driver's file; the name of the copy carries the set and the run,
because the three runs R1 (sets 2e and 2s, candidate terms at N2) have the same line in their header.

Check made before a file is written: d1 (N1 of the ligand to CZ of Phe131) and d2 (N1 to OD1 of Asn67) are
measured on the frame, with the atoms taken from `selections.json` of the built system (serial number = index
+ 1, as the header of the file says), and must equal the last line of the time series of the run
(md/results/series/, time 100.00 ns) within 0.002 A; the coordinates of a PDB file carry three decimals.

Written beside the frames: README.md, SOURCES.json and SHA256SUMS.  The run folders are read only.

    python3 md/production/make_final_frames.py [--runs DIR] [--out DIR] [--series DIR]
"""
import argparse
import csv
import datetime
import gzip
import hashlib
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_MD = os.path.dirname(HERE)
RUNS = [("bonded_2e", r) for r in ("X1", "X2", "X3", "N1", "N2", "N3", "R1", "R2", "R3")] + \
       [("bonded_2s", "R1"), ("bonded_2e_n2fit", "R1")]
TOL = 0.002                      # A: three decimals of the PDB coordinates against four of the series
N_ATOMS = 4094
IST = datetime.timezone(datetime.timedelta(hours=5, minutes=30))


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def sha256_file(path):
    with open(path, "rb") as fh:
        return sha256_bytes(fh.read())


def ist_now():
    t = datetime.datetime.now(IST)
    return "%d %s %d, %d:%02d %s IST" % (t.day, t.strftime("%b"), t.year, (t.hour % 12) or 12, t.minute,
                                         "AM" if t.hour < 12 else "PM")


def atoms_of(pdb_text):
    """{serial: (name, residue name, residue number, x, y, z)} of the ATOM lines."""
    out = {}
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            out[int(line[6:11])] = (line[12:16].strip(), line[17:20].strip(), int(line[22:26]),
                                    float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return out


def dist(a, b):
    return math.sqrt(sum((a[k] - b[k]) ** 2 for k in (3, 4, 5)))


def last_line_of_series(series_dir, model, run):
    with gzip.open(os.path.join(series_dir, "%s_%s.csv.gz" % (model, run)), "rt") as fh:
        rows = list(csv.DictReader(fh))
    return rows[-1]


def one_run(runs_dir, model, run, out_dir, series_dir):
    run_dir = os.path.join(runs_dir, model, run)
    src = os.path.join(run_dir, "final_solute.pdb")
    with open(src, "rb") as fh:
        raw = fh.read()
    text = raw.decode("ascii")
    with open(os.path.join(run_dir, "manifest.json")) as fh:
        man = json.load(fh)
    sel_path = os.path.join(man["system_dir"], "selections.json")
    with open(sel_path) as fh:
        sel = json.load(fh)
    at = atoms_of(text)
    if len(at) != N_ATOMS:
        raise SystemExit("%s/%s: %d atoms in final_solute.pdb, %d expected" % (model, run, len(at), N_ATOMS))
    n1, cz, od1 = (at[sel["lig_heavy"]["N1"] + 1], at[sel["phe131_CZ"] + 1], at[sel["asn67_OD1"] + 1])
    if (n1[0], cz[:2], od1[:2]) != ("N1", ("CZ", "PHE"), ("OD1", "ASN")):
        raise SystemExit("%s/%s: the atoms of selections.json are not N1, CZ of Phe and OD1 of Asn in the file"
                         % (model, run))
    d1, d2 = dist(n1, cz), dist(n1, od1)
    last = last_line_of_series(series_dir, model, run)
    if last["time_ns"] != "100.00":
        raise SystemExit("%s/%s: the series does not end at 100.00 ns" % (model, run))
    dd1, dd2 = d1 - float(last["d1_A"]), d2 - float(last["d2_A"])
    if abs(dd1) > TOL or abs(dd2) > TOL:
        raise SystemExit("%s/%s: the frame gives d1 %.4f and d2 %.4f, the series %s and %s at 100.00 ns"
                         % (model, run, d1, d2, last["d1_A"], last["d2_A"]))
    name = "%s_%s_final_solute.pdb.gz" % (model, run)
    out = os.path.join(out_dir, name)
    with open(out, "wb") as fh:
        with gzip.GzipFile(filename="", mode="wb", fileobj=fh, mtime=0, compresslevel=9) as gz:
            gz.write(raw)
    with gzip.open(out, "rb") as fh:
        if fh.read() != raw:
            raise SystemExit("%s: the compressed copy does not give back the source" % name)
    header = [ln for ln in text.splitlines() if ln.startswith(("REMARK", "CRYST1"))]
    return {"run": "%s/%s" % (model, run), "file": name, "bytes": os.path.getsize(out), "sha256": sha256_file(out),
            "source": src, "source_bytes": len(raw), "sha256_uncompressed_equal_to_the_source": sha256_bytes(raw),
            "n_atoms": len(at), "header_lines": header, "selections": sel_path,
            "selections_sha256": sha256_file(sel_path),
            "check_against_the_series": {"time_ns": last["time_ns"], "d1_frame_A": round(d1, 4),
                                         "d1_series_A": float(last["d1_A"]), "d2_frame_A": round(d2, 4),
                                         "d2_series_A": float(last["d2_A"]),
                                         "largest_difference_A": round(max(abs(dd1), abs(dd2)), 4),
                                         "state_of_the_series": last["state"]},
            "residue_numbers": {"first_and_last_deposited_CA": [sel["residue_numbers_deposited_for_CA"][0],
                                                                sel["residue_numbers_deposited_for_CA"][-1]],
                                "phe131_CZ_residue_in_file": cz[2], "asn67_OD1_residue_in_file": od1[2],
                                "ligand_residue_in_file": n1[2]}}


README = """# Last solute frame of the eleven runs of 100 ns

Written by `md/production/make_final_frames.py`. Each file is the `final_solute.pdb` that the driver
(`md/run_md.py`) wrote at the end of the last production segment of a run, compressed with gzip and otherwise
unchanged: the {n_atoms} atoms of the solute (protein with its acetyl cap, zinc, ligand) at 100 ns, coordinates in
Å as integrated, with the cell in the `CRYST1` line. Water and ions are not in these files. {n_files} files,
{total_kb} kB together. `SOURCES.json` names every source with its sha256, which is also the sha256 of the
uncompressed content of the copy; `SHA256SUMS` holds the sums of the files of this folder.

Check made for every file: d1 (N1 of the ligand to CZ of Phe131) and d2 (N1 to OD1 of Asn67) measured on the
frame equal the last line of the time series of the run (`md/results/series/`, 100.00 ns) within {worst:.4f} Å
(limit 0.002 Å; the coordinates of a PDB file carry three decimals).

## Numbering in the files

Atoms carry the serial number "index in the topology + 1". Residues are numbered in sequence: the acetyl cap is
residue 1 and His3 of the deposited structure 6RL9 is residue 2, so a residue with the deposited number *n*
is residue *n* − 1 up to residue 125 and *n* − 2 from residue 127 on (6RL9 has no residue 126; Phe131 is
residue {phe}, Asn67 residue {asn}). The zinc is residue 260 (`ZN1`) and the ligand residue {lig} (`SN1`; N1 is the
aniline nitrogen, N2 the sulfonamide nitrogen, `HN21` its hydrogen). Histidines are named by their tautomer
(`HID`, `HIE`); the three zinc ligands His94, His96 and His119 carry the names `HD1`, `HD2` and `HE1` of the
bonded zinc model.

## Files

| File | Run | Size (bytes) | d1; d2 on the frame (Å) | State at 100 ns |
|---|---|---|---|---|
{rows}

`bonded_2s_R1` and `bonded_2e_n2fit_R1` are the two sensitivity runs of start R (the other zinc set; the
candidate terms at N2, which gate G-R2 did not accept). They are in no count of rule 3 of the plan.
"""


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--runs", default="/root/cqc-sfa-md/runs")
    ap.add_argument("--out", default=os.path.join(REPO_MD, "results", "final_frames"))
    ap.add_argument("--series", default=os.path.join(REPO_MD, "results", "series"))
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    recs = [one_run(args.runs, m, r, args.out, args.series) for m, r in RUNS]
    worst = max(r["check_against_the_series"]["largest_difference_A"] for r in recs)
    rn = recs[0]["residue_numbers"]
    assert all(r["residue_numbers"] == rn for r in recs)
    rows = "\n".join("| `%s` | %s | %d | %.3f; %.3f | %s |"
                     % (r["file"], r["run"], r["bytes"], r["check_against_the_series"]["d1_frame_A"],
                        r["check_against_the_series"]["d2_frame_A"],
                        {"C": "crystal-like", "R": "rotated-like", "I": "intermediate"}[
                            r["check_against_the_series"]["state_of_the_series"]]) for r in recs)
    with open(os.path.join(args.out, "README.md"), "w", encoding="utf-8") as fh:
        fh.write(README.format(n_atoms="{:,}".format(N_ATOMS), n_files=len(recs),
                               total_kb="{:,.0f}".format(sum(r["bytes"] for r in recs) / 1000.0), worst=worst,
                               phe=rn["phe131_CZ_residue_in_file"], asn=rn["asn67_OD1_residue_in_file"],
                               lig=rn["ligand_residue_in_file"], rows=rows))
    with open(os.path.join(args.out, "SOURCES.json"), "w", encoding="utf-8") as fh:
        json.dump({"what": "Sources of the last solute frames of md/results/final_frames/",
                   "script": {"name": "md/production/make_final_frames.py",
                              "sha256": sha256_file(os.path.abspath(__file__))},
                   "written": ist_now(), "runs_dir": args.runs, "tolerance_A": TOL,
                   "largest_difference_from_the_series_A": worst, "runs": recs}, fh, indent=1)
        fh.write("\n")
    with open(os.path.join(args.out, "SHA256SUMS"), "w", encoding="ascii") as fh:
        for r in recs:
            fh.write("%s  %s\n" % (r["sha256"], r["file"]))
        for extra in ("README.md", "SOURCES.json"):
            fh.write("%s  %s\n" % (sha256_file(os.path.join(args.out, extra)), extra))
    print("written: %d frames in %s (%d bytes); largest difference from the series %.4f A"
          % (len(recs), args.out, sum(r["bytes"] for r in recs), worst))
    return 0


if __name__ == "__main__":
    sys.exit(main())
