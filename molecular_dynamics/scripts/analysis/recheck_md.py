#!/usr/bin/env python3
"""
recheck_md.py -- second, independent implementation of md/analysis/SPEC.md (version 1.3)

Stage S6 of md/PLAN.md: every number of the MD campaign is derived twice, by two programs written
from the specification by different authors, and the two must agree.  This is the second program.
It was written from md/analysis/SPEC.md alone; md/analysis/analyse_md.py and md/run_md.py were not
read before this program had produced its numbers for the test runs (version 1.2 on 6 Oct 2026,
9:04 AM IST; the changes of version 1.3, sections 9D and 9E and appendix A, on the same day after
11:30 AM IST, again from the text alone).  md/analysis/RECHECK.md says what was compared afterwards,
what was found, what was changed in this program after a comparison and the readings chosen where
the specification leaves a choice.

Version 1.3 of the specification: two water counts per full frame (section 9D).  The first is the
count that the driver acts on: the water oxygens within the cut-off c of the crystal positions of
the site atoms P, with P and c of the run (manifest hydration_check.settings when both are usable,
else C1 to C6 and 2.5 A).  The second is the wider count of version 1.2 (C1 to C6 and N1, 3.5 A).
The driver's record at the end of stage 03 and of every extension step is set beside the first.

Route, chosen to differ from the first program wherever the specification leaves the choice:
  - no MDAnalysis: the DCD files are read here (CHARMM layout, memory-mapped records, frames counted
    from the file size and not from the header), the PDB files by fixed columns;
  - the fit on the CA atoms uses Horn's quaternion method (largest eigenvector of a 4 x 4 matrix),
    not a singular value decomposition;
  - angles by 2 atan2(|u^ - v^|, |u^ + v^|) on unit vectors, torsions by the two plane normals
    (n1 = b1 x b2, n2 = b2 x b3, atan2(((n1 x n2) . b2/|b2|), n1 . n2)); both give the angles of
    section 0 of the specification;
  - minimum images through fractional coordinates of the cell of the frame (proved exact below half
    the smallest width of the cell; a search over 27 cells otherwise);
  - RMSF from the sums of the deviations from the reference positions.

Units: coordinates and distances in Angstrom, angles and torsions in degrees, times in ns, frame
intervals in ps.  No random numbers.  One thread (OMP_NUM_THREADS and friends are set to 1 below).

Use (PY = /root/cqc-sfa-tools/md/bin/python):

  PY recheck_md.py RUN_DIR [--system DIR] [--start X|N|R] [--model NAME] [--planned-ns NS]
                   [--solute-interval-ps PS] [--full-interval-ps PS] [--out DIR]
                   [--no-hash-trajectories] [--skip-name-check]
        analyse one run; writes timeseries_solute.csv, timeseries_full.csv, rmsf_ca.csv,
        summary.json (and the two equilibration series) to RUN_DIR/analysis_recheck/ (or --out)

  PY recheck_md.py --campaign DIR [--out DIR]
        analyse the production runs below DIR and evaluate the rules of the plan

  PY recheck_md.py --compare RUN_DIR [--primary DIR] [--out DIR] [options of a run]
        analyse RUN_DIR and compare with the output of analyse_md.py in RUN_DIR/analysis/ (or
        --primary) within the tolerances of section 14 of the specification

  PY recheck_md.py --compare-campaign DIR [--out DIR] [--run-primary]
        THE ONE COMMAND OF STAGE S6: every run folder below DIR (the validation run X0 included) is
        analysed here and compared with analyse_md.py's output; the campaign rules are compared with
        DIR/analysis_campaign/campaign.json.  --run-primary first calls analyse_md.py (as a separate
        process, it is never imported) where its output is missing or older than a trajectory file.
        Report: recheck_report.json and recheck_report.md in DIR/analysis_campaign_recheck/ (or --out).
        Exit status 0 = the two programs agree, 1 = at least one disagreement beyond tolerance,
        2 = something could not be compared (no output of the first program, or a run that this
        program cannot analyse), 3 = no disagreement, but a verdict, a part of gate G-V, a rule or a
        role differs between the two programs through frames at a threshold or a mean at a limit
        (section 14 allows that and asks for both readings: they are in the report and must be read).

What "agree" means (section 14 of the specification): every number within its tolerance, every count,
flag, state, labelling, sign, verdict, role and truth value identical, the same lines in the series,
the same keys in the compared blocks of summary.json, the same number of warnings.  A frame whose
value lies within the tolerance of a threshold may be classed differently by two correct programs;
such frames are listed with both readings and are not counted as disagreements.  What follows from
them is listed in the same way, but only where it can follow: a difference under a key of
summary.json is put down to threshold frames only if a column that the key depends on has such
frames in this run (THRESHOLD_SOURCES), a count may differ by at most their number, and the mean,
standard deviation and block means of a water count by at most what these frames can make.  A
verdict, a part of gate G-V, rule 1 or a role that differs in this way gives exit status 3.  Texts that the specification does not fix (notes, the names in
stats.<column>.kind, the label of a run) are listed where they differ and are not counted.

Written 6 Oct 2026.
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import csv
import hashlib
import json
import math
import platform
import re
import struct
import subprocess
import sys

import numpy as np

SPEC_VERSION = "1.3"
PROGRAM = "recheck_md.py"

# ----------------------------------------------------------------------------------------------
# constants of the specification (section in brackets)
# ----------------------------------------------------------------------------------------------
AKMA_PS = 0.04888821            # time unit of the step length in a DCD header, in ps
MID_D1_A = 6.7                  # [4]
MID_D2_A = 5.2                  # [4]
HB_DIST_A = 3.5                 # [5, 7]
HB_ANGLE_DEG = 135.0            # [5, 7]
FIFTH_A = 2.6                   # [6, 9C]
WATER_A = 3.5                   # [7] water at the aniline N; [9D] the cut-off of the wider count
RING_PLACE_ATOMS = ["C1", "C2", "C3", "C4", "C5", "C6"]         # [9D] defaults of the first count
RING_PLACE_CUTOFF_A = 2.5                                       # [9D]
RING_WIDE_ATOMS = ["C1", "C2", "C3", "C4", "C5", "C6", "N1"]    # [9D] the wider count: constants
SETTINGS_TOL_A = 1e-9           # [9D] |c - 2.5| for settings_are_the_defaults
STEP_END_TOL_PS = 1e-6          # [9E] |S + m dt - E| for "compared"
WINDOW_NS = 20.0                # [4] the last 20 ns: 2000 frames at 10 ps
DWELL_NS = 1.0                  # [4, 9A] 1 ns: 100 frames at 10 ps
BLOCK_NS = 10.0                 # [10] 10 ns: 1000 solute frames at 10 ps, 100 full frames at 100 ps
# The specification gives these three as frame counts "at" the default intervals and as times.  The
# plan (section 5, rule 2: "for the last 20 ns") fixes the time, so with another frame interval the
# frame counts follow from the time: the nearest whole number of frames, at least one.
RULE1_A = 2.0                   # [5]
VALIDATION_NS = 5.0             # [1]
DEFAULT_NS = 100.0
DEFAULT_SOLUTE_PS = 10.0
DEFAULT_FULL_PS = 100.0
REPEAT_A = 0.001                # [1] repeated frame
TWO_TRAJ_A = 0.01               # [1]
WHOLE_CA_A = 4.5                # [1]
GV1_MAX_A = 1.0                 # [11]
GV3A_WINDOW_A = (1.83, 2.03)
GV3B_WINDOW_A = (1.91, 2.12)
GV3C_TOL_DEG = 10.0
GV3C_REF_DEG = [("ang_N2_Zn_H94NE2_deg", 108.9), ("ang_N2_Zn_H96NE2_deg", 114.3),
                ("ang_N2_Zn_H119ND1_deg", 116.8), ("ang_H94NE2_Zn_H96NE2_deg", 103.8),
                ("ang_H94NE2_Zn_H119ND1_deg", 113.6), ("ang_H96NE2_Zn_H119ND1_deg", 98.4)]
GV4A_REF_DEG, GV4A_TOL_DEG = 108.8, 3.0
GV4B_REF_DEG, GV4B_TOL_DEG = 132.1, 20.0
EQ_CONSTANT_TABLE = [(1, 0.0, 100.0), (2, 100.0, 400.0), (3, 400.0, 700.0),
                     (4, 700.0, 1100.0), (5, 1100.0, 1500.0), (6, 1500.0, 2000.0)]   # [9E] ps
DRIVER_NEAREST_TOL_A = 0.002    # [9E]

LIG = ["C1", "C2", "C3", "C4", "C5", "C6", "N1", "N2", "O1", "O2", "S"]
LABELLINGS = ["identity", "ring_flip", "oxygen_swap", "ring_flip+oxygen_swap"]
SITE_NAMES = [("zn", "ZN"), ("his94_NE2", "NE2"), ("his96_NE2", "NE2"), ("his119_ND1", "ND1"),
              ("thr199_OG1", "OG1"), ("thr199_HG1", "HG1"), ("thr199_N", "N"), ("thr199_H", "H"),
              ("glu106_OE1", "OE1"), ("glu106_OE2", "OE2"), ("phe131_CZ", "CZ"),
              ("asn67_OD1", "OD1")]

SOLUTE_COLUMNS = [
    "lig_rmsd_A", "lig_rmsd_labelling", "lig_rmsd_identity_A", "d1_N1_F131CZ_A", "d2_N1_N67OD1_A",
    "state", "state_assigned", "d_N2_T199OG1_A", "ang_N2_HN2_T199OG1_deg", "hb_N2H_T199OG1",
    "d_Zn_N2_A", "d_Zn_H94NE2_A", "d_Zn_H96NE2_A", "d_Zn_H119ND1_A", "ang_N2_Zn_H94NE2_deg",
    "ang_N2_Zn_H96NE2_deg", "ang_N2_Zn_H119ND1_deg", "ang_H94NE2_Zn_H96NE2_deg",
    "ang_H94NE2_Zn_H119ND1_deg", "ang_H96NE2_Zn_H119ND1_deg", "ang_S_N2_Zn_deg", "d_Zn_O1_A",
    "d_Zn_O2_A", "d_Zn_Osulfonyl_min_A", "fifth_sulfonylO", "d_T199N_O1_A", "d_T199N_O2_A",
    "ang_T199N_H_O1_deg", "ang_T199N_H_O2_deg", "d_T199N_Osulfonyl_min_A", "hb_T199NH_Osulfonyl",
    "d_T199OG1_E106OE1_A", "d_T199OG1_E106OE2_A", "ang_T199OG1_HG1_E106OE1_deg",
    "ang_T199OG1_HG1_E106OE2_deg", "d_T199OG1_E106OE_min_A", "hb_T199OG1H_E106OE",
    "tors_C3_C4_S_N2_deg", "tors_C5_C4_S_N2_deg", "tau_fold_deg", "ang_C4_S_N2_deg", "d_S_N2_A",
    "ca_rmsd_A",
    "n2_volume", "n2_sign", "n2_sign_assigned", "d_Zn_E106OE1_A", "d_Zn_E106OE2_A",
    "d_Zn_E106OE_min_A", "fifth_E106OE1", "fifth_E106OE2", "fifth_E106OE"]
ND2_COLUMNS = ["d_N1_N67ND2_A", "asn67_ND2_nearer_N1"]
FULL_COLUMNS = [
    "d_Zn_Owater_min_A", "water_O_index_nearest_Zn", "d_Zn_Oprotein_min_A",
    "protein_O_index_nearest_Zn", "d_Zn_Osulfonyl_min_full_A", "d_Zn_O_any_min_A", "fifth_water",
    "fifth_proteinO", "fifth_any", "n_water_N1_3p5A", "d_Zn_N2_full_A", "n_water_ring_wide_3p5A",
    "d_water_ring_wide_min_A", "n_water_ring_place", "d_water_ring_place_min_A"]
FULL_COUNTS = ["n_water_N1_3p5A", "n_water_ring_wide_3p5A", "n_water_ring_place"]
SOLUTE_FLAGS = ["hb_N2H_T199OG1", "hb_T199NH_Osulfonyl", "hb_T199OG1H_E106OE", "fifth_sulfonylO",
                "fifth_E106OE1", "fifth_E106OE2", "fifth_E106OE"]
FULL_FLAGS = ["fifth_water", "fifth_proteinO", "fifth_any"]
FLAG_ORDER = ["hb_N2H_T199OG1", "hb_T199NH_Osulfonyl", "hb_T199OG1H_E106OE", "fifth_sulfonylO",
              "fifth_water", "fifth_proteinO", "fifth_any", "fifth_E106OE1", "fifth_E106OE2",
              "fifth_E106OE"]
TABLE_ROW_KEYS = [
    "run", "model", "start", "replica", "seen_ns", "planned_ns", "complete", "verdict",
    "provisional", "frac_start_state_last20ns", "frac_opposite_state_last20ns",
    "frac_intermediate_last20ns", "first_passage_ns", "n_transitions",
    "d1_N1_F131CZ_mean_last_half_A", "d2_N1_N67OD1_mean_last_half_A",
    "d_N1_N67ND2_mean_last_half_A", "asn67_ND2_nearer_N1_fraction", "lig_rmsd_mean_A",
    "lig_rmsd_sd_A", "lig_rmsd_mean_last_half_A", "lig_rmsd_sd_last_half_A", "rule1_left_crystal",
    "hb_N2H_T199OG1_fraction", "hb_T199NH_Osulfonyl_fraction", "hb_T199OG1H_E106OE_fraction",
    "d_Zn_N2_mean_A", "d_Zn_N2_sd_A", "d_Zn_H94NE2_mean_A", "d_Zn_H96NE2_mean_A",
    "d_Zn_H119ND1_mean_A", "ang_N_Zn_N_smallest_mean_deg", "ang_N_Zn_N_largest_mean_deg",
    "ang_N_Zn_N_largest_abs_deviation_deg", "d_Zn_Osulfonyl_min_mean_A", "fifth_sulfonylO_fraction",
    "d_Zn_Owater_min_mean_A", "d_Zn_Owater_min_min_A", "fifth_any_fraction", "n_water_N1_3p5A_mean",
    "tau_fold_mean_deg", "tau_fold_sd_deg", "ang_C4_S_N2_mean_deg", "ang_C4_S_N2_sd_deg",
    "d_S_N2_mean_A", "d_S_N2_sd_A", "ca_rmsd_mean_A", "ca_rmsd_mean_last_half_A", "rmsf_ca_max_A",
    "rmsf_ca_max_residue", "gate_GV_passed", "n_warnings",
    "n2_fraction_positive", "n2_fraction_positive_last20ns", "n2_n_lasting_sign_changes",
    "n2_first_lasting_sign_change_ns", "hb_N2H_T199OG1_fraction_before_first_passage",
    "hb_N2H_T199OG1_fraction_after_first_passage", "d_Zn_E106OE1_mean_A", "d_Zn_E106OE1_min_A",
    "fifth_E106OE1_fraction", "d_Zn_E106OE2_mean_A", "d_Zn_E106OE2_min_A", "fifth_E106OE2_fraction",
    "n_water_ring_wide_3p5A_mean", "d_water_ring_wide_min_min_A", "state_first_frame",
    "state_first_frame_differs_from_start", "n2_sign_first_frame", "equilibration_frames_solute",
    "n_water_ring_place_mean", "d_water_ring_place_min_min_A"]

SEG_RE = re.compile(r"^(solute|full)_seg(\d+)\.dcd$")
EQ_FILE_RE = re.compile(r"^(solute|full)_(eq\d+(?:x[1-9]\d*)?)\.dcd$")
EQ_TAG_RE = re.compile(r"^eq(\d+)(?:x([1-9]\d*))?$")
RUN_NAME_RE = re.compile(r"^([XNR])(\d+)$")
CAMPAIGN_NAME_RE = re.compile(r"^[XNR][1-9]\d*$")


class Stop(Exception):
    """The analysis cannot go on (a guard of section 1 of the specification that stops)."""


# ----------------------------------------------------------------------------------------------
# small helpers
# ----------------------------------------------------------------------------------------------
def r6(x):
    """A number as it is written to summary.json: six decimals; None for a missing value."""
    if x is None:
        return None
    x = float(x)
    if not math.isfinite(x):
        return None
    y = round(x, 6)
    return 0.0 if y == 0 else y


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for part in iter(lambda: fh.read(1 << 20), b""):
            h.update(part)
    return h.hexdigest()


def mget(d, *keys):
    for k in keys:
        if not isinstance(d, dict) or k not in d:
            return None
        d = d[k]
    return d


def is_number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def is_whole(x):
    return (isinstance(x, int) and not isinstance(x, bool)) or \
        (isinstance(x, float) and math.isfinite(x) and x == int(x))


def is_json_integer(x):
    """A JSON number written without fraction and exponent (Python reads it as int); not true/false."""
    return isinstance(x, int) and not isinstance(x, bool)


def ring_place_settings(manifest):
    """Section 9D, "Settings of the first count": the site atoms P and the cut-off c of the run.
    Both from manifest hydration_check.settings if both are usable there, else both from the
    defaults of the specification; never one of each.  Returns the block water_ring_place of
    summary.json and the warnings (at most one of each kind for a run)."""
    warn = []
    P, c, source = list(RING_PLACE_ATOMS), RING_PLACE_CUTOFF_A, "specification (defaults)"
    hc = manifest.get("hydration_check") if isinstance(manifest, dict) else None
    st = hc.get("settings") if isinstance(hc, dict) else None      # null counts as no key
    if st is not None:
        usable = False
        if isinstance(st, dict):
            names, cut = st.get("site_atoms"), st.get("water_cutoff_A")
            a_ok = (isinstance(names, list) and len(names) >= 1
                    and all(isinstance(x, str) and x in LIG for x in names)
                    and len(set(names)) == len(names))
            b_ok = (isinstance(cut, (int, float)) and not isinstance(cut, bool)
                    and math.isfinite(cut) and cut > 0)
            usable = a_ok and b_ok
        if usable:
            P, c, source = list(names), float(cut), "manifest hydration_check.settings"
        else:
            warn.append("hydration_check.settings of the manifest cannot be used (site_atoms must be a list of "
                        "different heavy-atom names of the ligand, water_cutoff_A a finite number above zero); "
                        "the defaults of the specification are used (C1 to C6, 2.5 A)")
    defaults = (len(P) == len(RING_PLACE_ATOMS) and set(P) == set(RING_PLACE_ATOMS)
                and abs(c - RING_PLACE_CUTOFF_A) <= SETTINGS_TOL_A)
    if source.startswith("manifest") and not defaults:
        warn.append("the run was made with a hydration check other than that of the plan (site atoms %s, "
                    "cut-off %g A); n_water_ring_place follows the run and is not the quantity of the plan's rule"
                    % (" ".join(P), c))
    block = {"note": "site atoms P and cut-off c of n_water_ring_place and d_water_ring_place_min_A "
                     "(section 9D); the wider count always uses C1 to C6 and N1 and 3.5 A",
             "site_atoms": P, "water_cutoff_A": c, "source": source, "settings_are_the_defaults": bool(defaults)}
    return block, warn


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return jsonable(o.tolist())
    if isinstance(o, np.str_):
        return str(o)
    return o


def write_json(path, obj):
    with open(path, "w") as fh:
        json.dump(jsonable(obj), fh, indent=1)
        fh.write("\n")


# ----------------------------------------------------------------------------------------------
# files: DCD and PDB
# ----------------------------------------------------------------------------------------------
class Dcd:
    """A DCD file (CHARMM layout).  Only whole frames count; their number comes from the file size."""

    def __init__(self, path):
        self.path = path
        self.name = os.path.basename(path)
        self.bytes = os.path.getsize(path)
        self.readable = False
        self.error = None
        self.n_frames = 0
        self.natoms = None
        self.first_step = None
        self.step_interval = None
        self.interval_ps = None
        self.first_frame_ps = None
        self.periodic = None
        self.trailing_bytes = 0
        self.header_n_frames = None
        self._mm = None
        try:
            self._header()
            self.readable = True
        except (ValueError, struct.error) as err:
            self.error = str(err)

    def _header(self):
        with open(self.path, "rb") as fh:
            b = fh.read(100)
            if len(b) < 100:
                raise ValueError("%d bytes: not readable as a DCD file" % self.bytes)
            en = None
            for e in ("<", ">"):
                if struct.unpack(e + "i", b[:4])[0] == 84:
                    en = e
                    break
            if en is None or b[4:8] != b"CORD":
                raise ValueError("no CORD record of 84 bytes: not a DCD file")
            ic = struct.unpack(en + "9i", b[8:44])
            delta = struct.unpack(en + "f", b[44:48])[0]
            fl = struct.unpack(en + "10i", b[48:88])
            if struct.unpack(en + "i", b[88:92])[0] != 84:
                raise ValueError("first record not closed: not a DCD file")
            tlen = struct.unpack(en + "i", b[92:96])[0]
            fh.seek(96 + tlen)
            t = fh.read(16)
            if len(t) < 16:
                raise ValueError("header cut off")
            tend, four, nat, four2 = struct.unpack(en + "4i", t)
            if tend != tlen or four != 4 or four2 != 4 or nat <= 0:
                raise ValueError("title or atom record damaged")
        if ic[8] != 0:
            raise ValueError("DCD file with fixed atoms: not supported")
        if fl[1] != 0:
            raise ValueError("DCD file with a fourth dimension: not supported")
        self.endian = en
        self.header_n_frames = ic[0]
        self.first_step = ic[1]
        self.step_interval = ic[2]
        self.delta_akma = float(delta)
        self.interval_ps = ic[2] * float(delta) * AKMA_PS
        self.first_frame_ps = ic[1] * float(delta) * AKMA_PS
        self.periodic = bool(fl[0])
        self.natoms = nat
        self.header_bytes = 96 + tlen + 16
        self.frame_bytes = (56 if self.periodic else 0) + 3 * (4 * nat + 8)
        body = self.bytes - self.header_bytes
        self.n_frames = body // self.frame_bytes
        self.trailing_bytes = body % self.frame_bytes
        if self.n_frames:
            mm = self._map()
            good = np.ones(self.n_frames, bool)
            for a in "xyz":
                good &= (mm[a + "0"] == 4 * nat) & (mm[a + "1"] == 4 * nat)
            if self.periodic:
                good &= (mm["c0"] == 48) & (mm["c1"] == 48)
            if not good.all():
                first_bad = int(np.argmin(good))
                self.damaged_from_frame = first_bad
                self.n_frames = first_bad
                self._mm = None

    def _map(self):
        if self._mm is None:
            en = self.endian
            f = []
            if self.periodic:
                f += [("c0", en + "i4"), ("cell", en + "f8", (6,)), ("c1", en + "i4")]
            for a in "xyz":
                f += [(a + "0", en + "i4"), (a, en + "f4", (self.natoms,)), (a + "1", en + "i4")]
            self._mm = np.memmap(self.path, dtype=np.dtype(f), mode="r", offset=self.header_bytes,
                                 shape=(self.n_frames,))
        return self._mm

    def read(self, i0, i1, idx=None):
        """Frames i0 .. i1-1: coordinates (F, n, 3) in double precision and the raw cell records
        (F, 6) in the order of the file (A, gamma, B, beta, alpha, C), or None."""
        mm = self._map()[i0:i1]
        cols = slice(None) if idx is None else np.asarray(idx)
        parts = [np.asarray(mm[a][:, cols], dtype=np.float64) for a in "xyz"]
        xyz = np.stack(parts, axis=-1)
        cell = np.array(mm["cell"], dtype=np.float64) if self.periodic else None
        return xyz, cell

    def info(self, index=None):
        d = {"file": self.name, "n_frames": int(self.n_frames), "bytes": int(self.bytes),
             "header_first_step": self.first_step, "header_step_interval": self.step_interval,
             "header_interval_ps": r6(self.interval_ps), "header_periodic": self.periodic,
             "readable": self.readable, "trailing_bytes": int(self.trailing_bytes)}
        if index is not None:
            d = dict([("index", index)] + list(d.items()))
        if self.error:
            d["error"] = self.error
        return d


def lattice(cell):
    """Raw cell records (F, 6) -> lattice vectors as rows H (F, 3, 3), their inverse, a flag per
    frame and half the smallest width of the cell.  An angle entry with |value| <= 1 is a cosine
    (CHARMM, OpenMM, the driver), otherwise degrees."""
    cell = np.atleast_2d(cell)
    A, g, B, b, a, C = [cell[:, k] for k in range(6)]
    ok = (A > 0) & (B > 0) & (C > 0) & np.isfinite(cell).all(axis=1)

    def cosine(v):
        return np.where(np.abs(v) <= 1.0, v, np.cos(np.radians(v)))
    ca, cb, cg = cosine(a), cosine(b), cosine(g)
    sg = np.sqrt(np.clip(1.0 - cg * cg, 1e-300, None))
    F = len(cell)
    H = np.zeros((F, 3, 3))
    H[:, 0, 0] = A
    H[:, 1, 0] = B * cg
    H[:, 1, 1] = B * sg
    H[:, 2, 0] = C * cb
    cy = (ca - cb * cg) / sg
    H[:, 2, 1] = C * cy
    H[:, 2, 2] = C * np.sqrt(np.clip(1.0 - cb * cb - cy * cy, 0.0, None))
    Hinv = np.zeros_like(H)
    half_width = np.zeros(F)
    for i in range(F):
        if ok[i] and abs(np.linalg.det(H[i])) > 1e-9:
            Hinv[i] = np.linalg.inv(H[i])
            half_width[i] = 0.5 / np.sqrt((Hinv[i] ** 2).sum(axis=0)).max()
        else:
            ok[i] = False
    return H, Hinv, ok, half_width


_OFFS = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)], float)
_ZERO = int(np.where((_OFFS == 0).all(axis=1))[0][0])


def nearest_lattice_vector(d, H, Hinv):
    """For difference vectors d (F, 3): the whole numbers n (F, 3) for which d - n H is shortest
    (nearest whole numbers of the fractional coordinates, then a search over the 27 neighbours)."""
    f = np.einsum("fi,fij->fj", d, Hinv)
    n0 = np.rint(f)
    cand = n0[:, None, :] + _OFFS[None, :, :]
    vec = d[:, None, :] - np.einsum("fki,fij->fkj", cand, H)
    d2 = (vec ** 2).sum(axis=-1)
    # prefer the rounded vector itself when it is as short as any other
    best = d2.argmin(axis=1)
    same = d2[np.arange(len(d)), _ZERO] <= d2.min(axis=1)
    best = np.where(same, _ZERO, best)
    return cand[np.arange(len(d)), best]


def image_distances(points, origin, H, Hinv, half_width, need_below):
    """Minimum-image distances of points (n, 3) from origin (3,) in one cell.  Rounding of the
    fractional coordinates is exact for every distance below half the smallest width of the cell;
    if the smallest distance found, or the largest distance the caller needs exactly, is not below
    that, all 27 neighbouring cells are searched."""
    d = points - origin
    if H is None:
        return np.sqrt((d * d).sum(axis=1))
    d = d - np.rint(d @ Hinv) @ H
    r2 = (d * d).sum(axis=1)
    if len(r2) and (math.sqrt(r2.min()) >= half_width or need_below >= half_width):
        best = r2.copy()
        for o in _OFFS:
            e = d - o @ H
            best = np.minimum(best, (e * e).sum(axis=1))
        r2 = best
    return np.sqrt(r2)


def read_pdb(path):
    """Atoms of ATOM and HETATM lines, in file order: names, residue names, residue numbers (text),
    coordinates."""
    names, resn, resi, xyz = [], [], [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                names.append(line[12:16].strip())
                resn.append(line[17:21].strip())
                resi.append(line[22:27].strip())
                try:
                    xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                except ValueError:
                    xyz.append((float("nan"),) * 3)
    return names, resn, resi, np.array(xyz, dtype=np.float64).reshape(-1, 3)


# ----------------------------------------------------------------------------------------------
# geometry
# ----------------------------------------------------------------------------------------------
def norm(v):
    return np.sqrt((v * v).sum(axis=-1))


def angle_deg(a, b, c):
    """Angle at b, 0 .. 180 degrees (Kahan's form; equal to atan2(|u x v|, u . v))."""
    u = a - b
    v = c - b
    u = u / norm(u)[..., None]
    v = v / norm(v)[..., None]
    return np.degrees(2.0 * np.arctan2(norm(u - v), norm(u + v)))


def cross(u, v):
    return np.stack([u[..., 1] * v[..., 2] - u[..., 2] * v[..., 1],
                     u[..., 2] * v[..., 0] - u[..., 0] * v[..., 2],
                     u[..., 0] * v[..., 1] - u[..., 1] * v[..., 0]], axis=-1)


def torsion_deg(a, b, c, d):
    """IUPAC torsion a-b-c-d in (-180, 180] degrees, from the normals of the two planes."""
    b1 = b - a
    b2 = c - b
    b3 = d - c
    n1 = cross(b1, b2)
    n2 = cross(b2, b3)
    y = (cross(n1, n2) * (b2 / norm(b2)[..., None])).sum(axis=-1)
    x = (n1 * n2).sum(axis=-1)
    t = np.degrees(np.arctan2(y, x))
    return np.where(t <= -180.0, 180.0, t)


def horn_fit(P, Q):
    """Least-squares superposition of P (F, M, 3) on Q (M, 3), equal weights, proper rotation.
    Returns R (F, 3, 3) and the centroids pbar (F, 3), qbar (3,):  x_fit = R (x - pbar) + qbar.
    Horn 1987: the unit quaternion is the eigenvector of the largest eigenvalue of N."""
    pbar = P.mean(axis=1)
    qbar = Q.mean(axis=0)
    Pc = P - pbar[:, None, :]
    Qc = Q - qbar
    S = np.einsum("fma,mb->fab", Pc, Qc)
    Sxx, Sxy, Sxz = S[:, 0, 0], S[:, 0, 1], S[:, 0, 2]
    Syx, Syy, Syz = S[:, 1, 0], S[:, 1, 1], S[:, 1, 2]
    Szx, Szy, Szz = S[:, 2, 0], S[:, 2, 1], S[:, 2, 2]
    N = np.empty((len(P), 4, 4))
    N[:, 0, 0] = Sxx + Syy + Szz
    N[:, 1, 1] = Sxx - Syy - Szz
    N[:, 2, 2] = -Sxx + Syy - Szz
    N[:, 3, 3] = -Sxx - Syy + Szz
    N[:, 0, 1] = N[:, 1, 0] = Syz - Szy
    N[:, 0, 2] = N[:, 2, 0] = Szx - Sxz
    N[:, 0, 3] = N[:, 3, 0] = Sxy - Syx
    N[:, 1, 2] = N[:, 2, 1] = Sxy + Syx
    N[:, 1, 3] = N[:, 3, 1] = Szx + Sxz
    N[:, 2, 3] = N[:, 3, 2] = Syz + Szy
    w, v = np.linalg.eigh(N)
    q = v[:, :, 3]
    q0, qx, qy, qz = q[:, 0], q[:, 1], q[:, 2], q[:, 3]
    R = np.empty((len(P), 3, 3))
    R[:, 0, 0] = q0 * q0 + qx * qx - qy * qy - qz * qz
    R[:, 0, 1] = 2 * (qx * qy - q0 * qz)
    R[:, 0, 2] = 2 * (qx * qz + q0 * qy)
    R[:, 1, 0] = 2 * (qy * qx + q0 * qz)
    R[:, 1, 1] = q0 * q0 - qx * qx + qy * qy - qz * qz
    R[:, 1, 2] = 2 * (qy * qz - q0 * qx)
    R[:, 2, 0] = 2 * (qz * qx - q0 * qy)
    R[:, 2, 1] = 2 * (qz * qy + q0 * qx)
    R[:, 2, 2] = q0 * q0 - qx * qx - qy * qy + qz * qz
    return R, pbar, qbar


# ----------------------------------------------------------------------------------------------
# statistics (section 10)
# ----------------------------------------------------------------------------------------------
def column_kind(col):
    if col.startswith("tors_"):
        return "tors"
    if col.endswith("_deg"):
        return "ang"
    if col.endswith("_A"):
        return "dist"
    if col == "n2_volume":
        return "vol"
    return "count"


def _mean_sd(v, circular):
    if len(v) == 0:
        return None, None
    if not circular:
        m = float(np.mean(v))
        return m, float(np.sqrt(np.mean((v - m) ** 2)))
    a = np.radians(v)
    s, c = float(np.sin(a).mean()), float(np.cos(a).mean())
    rr = math.hypot(s, c)
    if rr <= 0.0:
        return math.degrees(math.atan2(s, c)), None
    sd = math.degrees(math.sqrt(max(0.0, -2.0 * math.log(min(rr, 1.0)))))
    return math.degrees(math.atan2(s, c)), sd


def stat_block(v, block, trajectory, col):
    kind = column_kind(col)
    circ = kind == "tors"
    v = np.asarray(v, dtype=np.float64)
    n = len(v)
    h = n // 2
    mean, sd = _mean_sd(v, circ)
    mh, sh = _mean_sd(v[n - h:], circ) if h else (None, None)
    nb = n // block
    bm = [_mean_sd(v[k * block:(k + 1) * block], circ)[0] for k in range(nb)]
    sem = None
    if not circ and nb >= 2:
        sem = float(np.std(np.array(bm), ddof=1) / math.sqrt(nb))
    return {"kind": kind, "n": n, "mean": r6(mean), "sd": r6(sd),
            "min": r6(v.min()) if n else None, "max": r6(v.max()) if n else None,
            "n_last_half": h, "mean_last_half": r6(mh), "sd_last_half": r6(sh),
            "block_frames": block, "block_means": [r6(x) for x in bm], "block_sem": r6(sem),
            "trajectory": trajectory}


def fraction_block(flag, block):
    f = np.asarray(flag, dtype=np.int64)
    n = len(f)
    h = n // 2
    c = int(f.sum())
    ch = int(f[n - h:].sum()) if h else 0
    nb = n // block
    return {"n": n, "count": c, "fraction": r6(c / n) if n else None,
            "n_last_half": h, "count_last_half": ch, "fraction_last_half": r6(ch / h) if h else None,
            "block_frames": block,
            "block_fractions": [r6(int(f[k * block:(k + 1) * block].sum()) / block) for k in range(nb)]}


def carry(values, keep, start):
    """Assigned value: the value of the frame if it is one of `keep`, otherwise the one before it;
    before the first frame the value of the start."""
    out = []
    cur = start
    for x in values:
        if x in keep:
            cur = x
        out.append(cur)
    return out


def stretches(assigned):
    """Maximal runs of equal values: list of (value, first frame, number of frames)."""
    res = []
    s = 0
    n = len(assigned)
    for i in range(1, n + 1):
        if i == n or assigned[i] != assigned[s]:
            res.append((assigned[s], s, i - s))
            s = i
    return res


# ----------------------------------------------------------------------------------------------
# the system of a run: selections, reference, topologies
# ----------------------------------------------------------------------------------------------
class System:
    def __init__(self, system_dir):
        self.dir = system_dir
        self.sel_path = os.path.join(system_dir, "selections.json")
        self.ref_path = os.path.join(system_dir, "reference.pdb")
        for p in (self.sel_path, self.ref_path):
            if not os.path.isfile(p):
                raise Stop("system file missing: %s" % p)
        with open(self.sel_path) as fh:
            sel = json.load(fh)
        need = ["zn", "lig_heavy", "lig_HN2", "lig_all", "his94_NE2", "his96_NE2", "his119_ND1",
                "thr199_OG1", "thr199_HG1", "thr199_N", "thr199_H", "glu106_OE1", "glu106_OE2",
                "phe131_CZ", "asn67_OD1", "protein_CA", "protein_heavy", "water_O", "solute",
                "residue_numbers_deposited_for_CA"]
        missing = [k for k in need if k not in sel or sel[k] is None]
        if missing:
            raise Stop("selections.json lacks the keys %s" % ", ".join(missing))
        self.sel = sel
        self.atom = {k: int(sel[k]) for k, _ in SITE_NAMES}
        for nm in LIG:
            if nm not in sel["lig_heavy"]:
                raise Stop("selections.json: lig_heavy lacks %s" % nm)
            self.atom[nm] = int(sel["lig_heavy"][nm])
        self.atom["HN2"] = int(sel["lig_HN2"])
        self.has_nd2 = sel.get("asn67_ND2") is not None
        if self.has_nd2:
            self.atom["asn67_ND2"] = int(sel["asn67_ND2"])
        self.lig_all = [int(i) for i in sel["lig_all"]]
        self.ca = [int(i) for i in sel["protein_CA"]]
        self.protein_heavy = [int(i) for i in sel["protein_heavy"]]
        self.water_O = [int(i) for i in sel["water_O"]]
        self.solute = [int(i) for i in sel["solute"]]
        self.resnum_ca = [int(i) for i in sel["residue_numbers_deposited_for_CA"]]
        if len(self.resnum_ca) != len(self.ca):
            raise Stop("selections.json: protein_CA and residue_numbers_deposited_for_CA differ in length")
        # reference
        names, resn, resi, xyz = read_pdb(self.ref_path)
        ca_rows = [i for i, nm in enumerate(names) if nm == "CA"]
        if len(ca_rows) != len(self.ca):
            raise Stop("reference.pdb has %d CA atoms, selections.json %d" % (len(ca_rows), len(self.ca)))
        self.ref_ca = xyz[ca_rows]
        self.ref_ca_resn = [resn[i] for i in ca_rows]
        ref_lig = []
        for nm in LIG:
            rows = [i for i, x in enumerate(names) if x == nm]
            if len(rows) != 1:
                raise Stop("reference.pdb has %d atoms named %s (one expected)" % (len(rows), nm))
            ref_lig.append(xyz[rows[0]])
        self.ref_lig = np.array(ref_lig)
        self.ref_by_name = {nm: self.ref_lig[i] for i, nm in enumerate(LIG)}
        # the four labellings: ref position that atom a is compared with
        k = {nm: i for i, nm in enumerate(LIG)}

        def perm(pairs):
            p = list(range(len(LIG)))
            for a, b in pairs:
                p[k[a]], p[k[b]] = k[b], k[a]
            return p
        self.perms = [perm([]), perm([("C2", "C6"), ("C3", "C5")]), perm([("O1", "O2")]),
                      perm([("C2", "C6"), ("C3", "C5"), ("O1", "O2")])]
        # neighbouring CA atoms: numbers i, i+1 and 125-127
        self.ca_pairs = [i for i in range(len(self.ca) - 1)
                         if self.resnum_ca[i + 1] - self.resnum_ca[i] == 1
                         or (self.resnum_ca[i], self.resnum_ca[i + 1]) == (125, 127)]

    def name_mismatches(self, names, resn, where):
        """Guard 'atom names'.  `where` maps an atom index of the topology to its position in the
        PDB file (None if the atom is not in that file)."""
        bad = []

        def check(label, index, test, want):
            p = where(index)
            if p is None:
                return
            if p >= len(names) or not test(names[p]):
                bad.append("%s (index %d): %s, expected %s" %
                           (label, index, names[p] if p < len(names) else "no atom", want))
        for key, want in SITE_NAMES:
            check(key, self.atom[key], lambda x, w=want: x == w, want)
        for nm in LIG:
            check("lig_heavy." + nm, self.atom[nm], lambda x, w=nm: x == w, nm)
        check("lig_HN2", self.atom["HN2"], lambda x: x.startswith("H"), "a name starting with H")
        if self.has_nd2:
            check("asn67_ND2", self.atom["asn67_ND2"], lambda x: x == "ND2", "ND2")
        for j, i in enumerate(self.ca):
            p = where(i)
            if p is None:
                continue
            check("protein_CA[%d]" % j, i, lambda x: x == "CA", "CA")
            if p < len(names) and resn[p] != self.ref_ca_resn[j]:
                bad.append("protein_CA[%d] (index %d): residue %s, reference.pdb has %s" %
                           (j, i, resn[p], self.ref_ca_resn[j]))
        n_bad_water = 0
        for i in self.water_O:
            p = where(i)
            if p is None:
                continue
            if p >= len(names) or not names[p].startswith("O"):
                n_bad_water += 1
        if n_bad_water:
            bad.append("%d atoms of water_O do not have a name starting with O" % n_bad_water)
        return bad


# ----------------------------------------------------------------------------------------------
# solute trajectory: every quantity of a frame (sections 2 to 9C)
# ----------------------------------------------------------------------------------------------
class SoluteReader:
    def __init__(self, system):
        s = self.system = system
        pos = {a: i for i, a in enumerate(s.solute)}
        needed = [s.atom[k] for k, _ in SITE_NAMES] + [s.atom[n] for n in LIG] + [s.atom["HN2"]] \
            + s.lig_all + s.ca + ([s.atom["asn67_ND2"]] if s.has_nd2 else [])
        absent = sorted({a for a in needed if a not in pos})
        if absent:
            raise Stop("atoms needed by the analysis are not in the list 'solute': %s" % absent[:10])
        cols = sorted({pos[a] for a in needed})
        self.read_cols = np.array(cols)
        sub = {p: k for k, p in enumerate(cols)}
        self.c = {key: sub[pos[idx]] for key, idx in s.atom.items()}
        self.c_lig_all = np.array([sub[pos[a]] for a in s.lig_all])
        self.c_ca = np.array([sub[pos[a]] for a in s.ca])
        self.c_lig = np.array([self.c[n] for n in LIG])
        self.pos = pos

    def new_accumulator(self):
        return {"cols": {}, "frames_with_box": 0, "frames_without_box": 0, "frames_zinc_moved": 0,
                "frames_ligand_moved": 0, "largest_neighbour_CA_distance_A": 0.0,
                "sum_dev": np.zeros((len(self.system.ca), 3)),
                "sum_dev2": np.zeros(len(self.system.ca)), "n": 0,
                "labelling_ties": [], "file_of_frame": []}

    def read_file(self, dcd, acc, tag, keep_raw_ca=False, chunk=500):
        """Add all frames of one DCD file to the accumulator.  Returns the raw CA coordinates of
        its frames (single precision) if asked."""
        raw = []
        for i0 in range(0, dcd.n_frames, chunk):
            i1 = min(dcd.n_frames, i0 + chunk)
            X, cell = dcd.read(i0, i1, self.read_cols)
            if keep_raw_ca:
                raw.append(X[:, self.c_ca].astype(np.float32))
            out = self.quantities(X, cell, acc)
            for k, v in out.items():
                acc["cols"].setdefault(k, []).append(v)
            acc["file_of_frame"].extend([tag] * (i1 - i0))
        return np.concatenate(raw) if raw else None

    def quantities(self, X, cell, acc):
        s = self.system
        c = self.c
        F = len(X)
        # --- section 2, step 1: images -------------------------------------------------------
        if cell is not None:
            H, Hinv, ok, _ = lattice(cell)
        else:
            ok = np.zeros(F, bool)
        nbox = int(ok.sum())
        acc["frames_with_box"] += nbox
        acc["frames_without_box"] += F - nbox
        if nbox:
            w = np.where(ok)[0]
            n = nearest_lattice_vector(X[w, c["zn"]] - X[w, c["his94_NE2"]], H[w], Hinv[w])
            X[w, c["zn"]] -= np.einsum("fi,fij->fj", n, H[w])
            acc["frames_zinc_moved"] += int((n != 0).any(axis=1).sum())
            n = nearest_lattice_vector(X[w, c["N2"]] - X[w, c["zn"]], H[w], Hinv[w])
            shift = np.einsum("fi,fij->fj", n, H[w])
            lig_rows = X[w][:, self.c_lig_all] - shift[:, None, :]
            X[np.ix_(w, self.c_lig_all)] = lig_rows
            acc["frames_ligand_moved"] += int((n != 0).any(axis=1).sum())

        def g(k):
            return X[:, c[k]]
        o = {}
        # --- section 2, step 2: fit on the CA atoms ------------------------------------------
        P = X[:, self.c_ca]
        R, pbar, qbar = horn_fit(P, s.ref_ca)
        ca_fit = np.einsum("fij,fmj->fmi", R, P - pbar[:, None, :]) + qbar
        dev = ca_fit - s.ref_ca
        o["ca_rmsd_A"] = np.sqrt((dev * dev).sum(axis=-1).mean(axis=-1))
        acc["sum_dev"] += dev.sum(axis=0)
        acc["sum_dev2"] += (dev * dev).sum(axis=-1).sum(axis=0)
        acc["n"] += F
        if s.ca_pairs:
            pr = np.array(s.ca_pairs)
            dn = norm(P[:, pr + 1] - P[:, pr]).max()
            acc["largest_neighbour_CA_distance_A"] = max(acc["largest_neighbour_CA_distance_A"], float(dn))
        # --- section 3: ligand RMSD, four labellings -----------------------------------------
        lig_fit = np.einsum("fij,fmj->fmi", R, X[:, self.c_lig] - pbar[:, None, :]) + qbar
        msd = np.stack([((lig_fit - s.ref_lig[p]) ** 2).sum(axis=-1).mean(axis=-1) for p in s.perms], axis=1)
        rm = np.sqrt(msd)
        lab = msd.argmin(axis=1)
        o["lig_rmsd_A"] = rm[np.arange(F), lab]
        o["lig_rmsd_labelling"] = lab
        o["lig_rmsd_identity_A"] = rm[:, 0]
        srt = np.sort(rm, axis=1)
        o["_labelling_gap_A"] = srt[:, 1] - srt[:, 0]
        # --- section 4: orientation ----------------------------------------------------------
        d1 = norm(g("N1") - g("phe131_CZ"))
        d2 = norm(g("N1") - g("asn67_OD1"))
        o["d1_N1_F131CZ_A"] = d1
        o["d2_N1_N67OD1_A"] = d2
        st = np.full(F, "I", dtype="<U1")
        st[(d1 < MID_D1_A) & (d2 > MID_D2_A)] = "C"
        st[(d1 > MID_D1_A) & (d2 < MID_D2_A)] = "R"
        o["state"] = st
        # --- section 5: N2-H...OG1(Thr199) ---------------------------------------------------
        d = norm(g("N2") - g("thr199_OG1"))
        a = angle_deg(g("N2"), g("HN2"), g("thr199_OG1"))
        o["d_N2_T199OG1_A"] = d
        o["ang_N2_HN2_T199OG1_deg"] = a
        o["hb_N2H_T199OG1"] = ((d <= HB_DIST_A) & (a >= HB_ANGLE_DEG)).astype(int)
        # --- section 6: zinc site ------------------------------------------------------------
        four = [("N2", "N2"), ("H94NE2", "his94_NE2"), ("H96NE2", "his96_NE2"), ("H119ND1", "his119_ND1")]
        for nm, key in four:
            o["d_Zn_%s_A" % nm] = norm(g("zn") - g(key))
        for i in range(4):
            for j in range(i + 1, 4):
                o["ang_%s_Zn_%s_deg" % (four[i][0], four[j][0])] = angle_deg(g(four[i][1]), g("zn"), g(four[j][1]))
        o["ang_S_N2_Zn_deg"] = angle_deg(g("S"), g("N2"), g("zn"))
        o["d_Zn_O1_A"] = norm(g("zn") - g("O1"))
        o["d_Zn_O2_A"] = norm(g("zn") - g("O2"))
        o["d_Zn_Osulfonyl_min_A"] = np.minimum(o["d_Zn_O1_A"], o["d_Zn_O2_A"])
        o["fifth_sulfonylO"] = (o["d_Zn_Osulfonyl_min_A"] <= FIFTH_A).astype(int)
        # --- section 7: other contacts -------------------------------------------------------
        da = norm(g("thr199_N") - g("O1"))
        db = norm(g("thr199_N") - g("O2"))
        aa = angle_deg(g("thr199_N"), g("thr199_H"), g("O1"))
        ab = angle_deg(g("thr199_N"), g("thr199_H"), g("O2"))
        o["d_T199N_O1_A"], o["d_T199N_O2_A"] = da, db
        o["ang_T199N_H_O1_deg"], o["ang_T199N_H_O2_deg"] = aa, ab
        o["d_T199N_Osulfonyl_min_A"] = np.minimum(da, db)
        o["hb_T199NH_Osulfonyl"] = (((da <= HB_DIST_A) & (aa >= HB_ANGLE_DEG)) |
                                    ((db <= HB_DIST_A) & (ab >= HB_ANGLE_DEG))).astype(int)
        da = norm(g("thr199_OG1") - g("glu106_OE1"))
        db = norm(g("thr199_OG1") - g("glu106_OE2"))
        aa = angle_deg(g("thr199_OG1"), g("thr199_HG1"), g("glu106_OE1"))
        ab = angle_deg(g("thr199_OG1"), g("thr199_HG1"), g("glu106_OE2"))
        o["d_T199OG1_E106OE1_A"], o["d_T199OG1_E106OE2_A"] = da, db
        o["ang_T199OG1_HG1_E106OE1_deg"], o["ang_T199OG1_HG1_E106OE2_deg"] = aa, ab
        o["d_T199OG1_E106OE_min_A"] = np.minimum(da, db)
        o["hb_T199OG1H_E106OE"] = (((da <= HB_DIST_A) & (aa >= HB_ANGLE_DEG)) |
                                   ((db <= HB_DIST_A) & (ab >= HB_ANGLE_DEG))).astype(int)
        # --- section 8: ligand internal coordinates ------------------------------------------
        t3 = torsion_deg(g("C3"), g("C4"), g("S"), g("N2"))
        t5 = torsion_deg(g("C5"), g("C4"), g("S"), g("N2"))
        o["tors_C3_C4_S_N2_deg"], o["tors_C5_C4_S_N2_deg"] = t3, t5
        o["tau_fold_deg"] = np.maximum(np.abs(t3), np.abs(t5))
        o["ang_C4_S_N2_deg"] = angle_deg(g("C4"), g("S"), g("N2"))
        o["d_S_N2_A"] = norm(g("S") - g("N2"))
        # --- section 9A: signed volume at N2 -------------------------------------------------
        ez = g("zn") - g("N2")
        eh = g("HN2") - g("N2")
        es = g("S") - g("N2")
        ez = ez / norm(ez)[:, None]
        eh = eh / norm(eh)[:, None]
        es = es / norm(es)[:, None]
        vol = np.linalg.det(np.stack([ez, eh, es], axis=1))
        o["n2_volume"] = vol
        o["n2_sign"] = np.where(vol > 0, 1, np.where(vol < 0, -1, 0)).astype(int)
        # --- section 9C: zinc to Glu106 ------------------------------------------------------
        e1 = norm(g("zn") - g("glu106_OE1"))
        e2 = norm(g("zn") - g("glu106_OE2"))
        o["d_Zn_E106OE1_A"], o["d_Zn_E106OE2_A"] = e1, e2
        o["d_Zn_E106OE_min_A"] = np.minimum(e1, e2)
        o["fifth_E106OE1"] = (e1 <= FIFTH_A).astype(int)
        o["fifth_E106OE2"] = (e2 <= FIFTH_A).astype(int)
        o["fifth_E106OE"] = (o["d_Zn_E106OE_min_A"] <= FIFTH_A).astype(int)
        # --- section 4, optional: Asn67 ND2 --------------------------------------------------
        if s.has_nd2:
            dn2 = norm(g("N1") - g("asn67_ND2"))
            o["d_N1_N67ND2_A"] = dn2
            o["asn67_ND2_nearer_N1"] = (dn2 < d2).astype(int)
        return o

    @staticmethod
    def finish(acc):
        cols = {}
        for k, parts in acc["cols"].items():
            cols[k] = np.concatenate(parts) if parts else np.array([])
        return cols


# ----------------------------------------------------------------------------------------------
# full trajectory: sections 6, 7 and 9D
# ----------------------------------------------------------------------------------------------
class FullReader:
    def __init__(self, system, names, place_atoms, place_cutoff):
        s = self.system = system
        self.water = np.array(s.water_O, dtype=np.int64)
        self.prot_O = np.array([i for i in s.protein_heavy if i < len(names) and names[i].startswith("O")],
                               dtype=np.int64)
        self.ca = np.array(s.ca, dtype=np.int64)
        # section 9D: the crystal positions needed, once each; the two counts are minima over subsets
        self.place_atoms = list(place_atoms)
        self.place_cutoff = float(place_cutoff)
        self.site_names = list(RING_WIDE_ATOMS) + [a for a in self.place_atoms if a not in RING_WIDE_ATOMS]
        self.site_ref = np.array([s.ref_by_name[a] for a in self.site_names])
        self.rows_wide = [self.site_names.index(a) for a in RING_WIDE_ATOMS]
        self.rows_place = [self.site_names.index(a) for a in self.place_atoms]

    def read_file(self, dcd, rows, aux, tag):
        s = self.system
        a = s.atom
        for i in range(dcd.n_frames):
            X, cell = dcd.read(i, i + 1)
            X = X[0]
            H = Hinv = None
            hw = float("inf")
            if cell is not None:
                Hs, Hinvs, ok, hws = lattice(cell)
                if ok[0]:
                    H, Hinv, hw = Hs[0], Hinvs[0], float(hws[0])
            if H is None:
                aux["frames_without_box"] += 1
            zn = X[a["zn"]]
            row = {"_tag": tag}
            # nearest water oxygen to the zinc
            if len(self.water):
                dw = image_distances(X[self.water], zn, H, Hinv, hw, FIFTH_A)
                k = int(dw.argmin())
                row["d_Zn_Owater_min_A"] = float(dw[k])
                row["water_O_index_nearest_Zn"] = int(self.water[k])
            else:
                row["d_Zn_Owater_min_A"] = float("nan")
                row["water_O_index_nearest_Zn"] = -1
            if len(self.prot_O):
                dp = image_distances(X[self.prot_O], zn, H, Hinv, hw, FIFTH_A)
                k = int(dp.argmin())
                row["d_Zn_Oprotein_min_A"] = float(dp[k])
                row["protein_O_index_nearest_Zn"] = int(self.prot_O[k])
            else:
                row["d_Zn_Oprotein_min_A"] = float("nan")
                row["protein_O_index_nearest_Zn"] = -1
            ds = image_distances(X[[a["O1"], a["O2"]]], zn, H, Hinv, hw, FIFTH_A)
            row["d_Zn_Osulfonyl_min_full_A"] = float(ds.min())
            three = [row["d_Zn_Owater_min_A"], row["d_Zn_Oprotein_min_A"], row["d_Zn_Osulfonyl_min_full_A"]]
            row["d_Zn_O_any_min_A"] = float(np.nanmin(three))
            row["fifth_water"] = int(row["d_Zn_Owater_min_A"] <= FIFTH_A)
            row["fifth_proteinO"] = int(row["d_Zn_Oprotein_min_A"] <= FIFTH_A)
            row["fifth_any"] = int(row["d_Zn_O_any_min_A"] <= FIFTH_A)
            # water at the aniline N
            if len(self.water):
                dn = image_distances(X[self.water], X[a["N1"]], H, Hinv, hw, WATER_A + 0.01)
                row["n_water_N1_3p5A"] = int((dn <= WATER_A).sum())
                edge_n1 = int((np.abs(dn - WATER_A) <= 0.001).sum())
            else:
                row["n_water_N1_3p5A"] = 0
                edge_n1 = 0
            row["d_Zn_N2_full_A"] = float(image_distances(X[[a["N2"]]], zn, H, Hinv, hw, 0.0)[0])
            # section 9D: water in the place of the crystal ring.  The CA atoms of the full frame are
            # fitted on those of the reference, x_fit = R (x - pbar) + qbar; the crystal positions are
            # carried back into the frame, s_a = R^T (r_a - qbar) + pbar (no atom is moved before).
            R, pbar, qbar = horn_fit(X[self.ca][None, :, :], s.ref_ca)
            sites = (self.site_ref - qbar) @ R[0] + pbar[0]
            c = self.place_cutoff
            if len(self.water):
                W = X[self.water]
                need = max(WATER_A, c) + 0.01
                per_site = np.array([image_distances(W, p, H, Hinv, hw, need) for p in sites])
                d_wide = per_site[self.rows_wide].min(axis=0)       # d_w of the wider count
                d_place = per_site[self.rows_place].min(axis=0)     # d_w of the count of the driver
                row["n_water_ring_wide_3p5A"] = int((d_wide <= WATER_A).sum())
                row["d_water_ring_wide_min_A"] = float(d_wide.min())
                row["n_water_ring_place"] = int((d_place <= c).sum())
                row["d_water_ring_place_min_A"] = float(d_place.min())
                edge_wide = int((np.abs(d_wide - WATER_A) <= 0.001).sum())
                edge_place = int((np.abs(d_place - c) <= 0.001).sum())
            else:
                row["n_water_ring_wide_3p5A"] = 0
                row["d_water_ring_wide_min_A"] = float("nan")
                row["n_water_ring_place"] = 0
                row["d_water_ring_place_min_A"] = float("nan")
                edge_wide = edge_place = 0
            row["_edge_n1"] = edge_n1
            row["_edge_wide"] = edge_wide
            row["_edge_place"] = edge_place
            rows.append(row)


# ----------------------------------------------------------------------------------------------
# who the run is (section 1)
# ----------------------------------------------------------------------------------------------
def identify_run(run_dir, args, manifest):
    folder = os.path.basename(os.path.normpath(run_dir))
    parent = os.path.basename(os.path.dirname(os.path.normpath(run_dir)))
    src = {}
    fm = RUN_NAME_RE.match(folder)
    rid = manifest.get("run_id") if isinstance(manifest, dict) else None
    rm = RUN_NAME_RE.match(rid) if isinstance(rid, str) else None
    # start
    start = None
    for label, val in (("command line (--start)", getattr(args, "start", None)),
                       ("manifest start", mget(manifest, "start")),
                       ("manifest build_json.start", mget(manifest, "build_json", "start")),
                       ("manifest run_id", rm.group(1) if rm else None),
                       ("folder name", fm.group(1) if fm else None)):
        if val in ("X", "N", "R"):
            start, src["start"] = val, label
            break
    if start is None:
        raise Stop("the start of the run cannot be found (no --start, no manifest entry, and the folder "
                   "name '%s' is not X, N or R followed by a number)" % folder)
    # replica
    replica = None
    src["replica"] = "none"
    mrep = mget(manifest, "replica")
    if mrep is not None:
        replica, src["replica"] = mrep, "manifest replica"
    elif rm:
        replica, src["replica"] = int(rm.group(2)), "manifest run_id"
    elif fm:
        replica, src["replica"] = int(fm.group(2)), "folder name"
    replica_is_int = isinstance(replica, int) and not isinstance(replica, bool)
    # model
    model = None
    for label, val in (("command line (--model)", getattr(args, "model", None)),
                       ("manifest model", mget(manifest, "model")),
                       ("manifest build_json.model", mget(manifest, "build_json", "model")),
                       ("parent folder", parent)):
        if isinstance(val, str) and val:
            model, src["model"] = val, label
            break
    role = "validation" if (start == "X" and replica_is_int and replica == 0) else "production"
    # planned length
    warnings = []
    planned = None
    cli_planned = getattr(args, "planned_ns", None)
    if cli_planned is not None:
        planned, src["planned_ns"] = float(cli_planned), "command line (--planned-ns)"
    elif role == "validation":
        planned, src["planned_ns"] = VALIDATION_NS, "specification (validation run: 5 ns)"
        named = mget(manifest, "planned_ns")
        if not is_number(named):
            named = mget(manifest, "production", "target_ns")
        if is_number(named) and abs(named - VALIDATION_NS) > 1e-9:
            warnings.append("validation run: the manifest names a length of %g ns; the planned length "
                            "stays 5 ns (specification, section 1)" % named)
    else:
        for label, val in (("manifest planned_ns", mget(manifest, "planned_ns")),
                           ("manifest production.target_ns", mget(manifest, "production", "target_ns")),
                           ("manifest protocol.production.length_ns",
                            mget(manifest, "protocol", "production", "length_ns"))):
            if is_number(val):
                planned, src["planned_ns"] = float(val), label
                break
        if planned is None:
            planned, src["planned_ns"] = DEFAULT_NS, "default (100 ns)"
    # intervals
    def interval(cli, plain, nested, default, key):
        for label, val in (("command line", cli), ("manifest " + plain, mget(manifest, plain)),
                           ("manifest protocol.output." + nested, mget(manifest, "protocol", "output", nested))):
            if is_number(val) and val > 0:
                src[key] = label
                return float(val)
        src[key] = "default (%g ps)" % default
        return default
    dt_s = interval(getattr(args, "solute_interval_ps", None), "solute_interval_ps",
                    "solute_frame_interval_ps", DEFAULT_SOLUTE_PS, "solute_interval_ps")
    dt_f = interval(getattr(args, "full_interval_ps", None), "full_interval_ps",
                    "full_frame_interval_ps", DEFAULT_FULL_PS, "full_interval_ps")
    # warnings of the role and of the names
    if fm:
        role_by_name = "validation" if (fm.group(1) == "X" and int(fm.group(2)) == 0) else "production"
        if role_by_name != role:
            warnings.append("the folder name %s alone would give the role %s; the sources give %s, which stands"
                            % (folder, role_by_name, role))
        if fm.group(1) != start:
            warnings.append("the start %s (%s) is not the letter of the folder name %s"
                            % (start, src["start"], folder))
    # system folder
    sysdir = None
    if getattr(args, "system", None):
        sysdir, src["system_dir"] = args.system, "command line (--system)"
    else:
        msys = mget(manifest, "system_dir")
        if isinstance(msys, str) and msys:
            cand = msys if os.path.isabs(msys) else os.path.join(run_dir, msys)
            if os.path.isdir(cand):
                sysdir, src["system_dir"] = cand, "manifest system_dir"
        if sysdir is None:
            sysdir = os.path.join(run_dir, "..", "..", "..", "systems", model, start)
            src["system_dir"] = "convention" + (" (the system_dir of the manifest is not on this machine)"
                                                if isinstance(msys, str) and msys else "")
    sysdir = os.path.normpath(sysdir)
    if replica_is_int:
        label = "%s/%s%d" % (model, start, replica)
    else:
        label = "%s/%s" % (model, folder)
    return {"folder": folder, "start": start, "replica": replica, "replica_is_int": replica_is_int,
            "model": model, "role": role, "planned_ns": planned, "dt_s": dt_s, "dt_f": dt_f,
            "sources": src, "system_dir": sysdir, "label": label, "warnings": warnings,
            "cli_planned": cli_planned is not None}


def header_check(dcds, dt_ps):
    """Guards 'interval' and 'contiguous headers' of section 1 for the segments of one kind."""
    with_frames = [d for d in dcds if d.readable and d.n_frames > 0]
    details = []
    interval = None
    if with_frames:
        interval = all(abs(d.interval_ps - dt_ps) <= 0.001 * dt_ps for d in with_frames)
        for d in with_frames:
            if abs(d.interval_ps - dt_ps) > 0.001 * dt_ps:
                details.append("%s: header interval %.6g ps, expected %.6g ps" % (d.name, d.interval_ps, dt_ps))
    contiguous = None
    if len(with_frames) >= 2 and not all(d.first_step == d.step_interval for d in with_frames):
        contiguous = True
        for a, b in zip(with_frames[:-1], with_frames[1:]):
            t_last = (a.first_step + (a.n_frames - 1) * a.step_interval) * a.delta_akma * AKMA_PS
            t_next = b.first_frame_ps
            if abs(t_next - (t_last + a.interval_ps)) > 0.001 * dt_ps:
                contiguous = False
                details.append("%s starts at %.6g ps by its header; %s ends at %.6g ps and one interval later is %.6g ps"
                               % (b.name, t_next, a.name, t_last, t_last + a.interval_ps))
    return {"interval_matches": interval, "contiguous": contiguous, "details": details}


# ----------------------------------------------------------------------------------------------
# one run
# ----------------------------------------------------------------------------------------------
def fmt_value(col, v):
    if isinstance(v, str):
        return v
    if col.endswith("_A") or col == "n2_volume":
        return "%.4f" % v
    if col.endswith("_deg"):
        return "%.3f" % v
    return "%d" % int(v)


def write_solute_csv(path, cols, head, lead_rows, columns):
    with open(path, "w") as fh:
        fh.write(",".join(head + columns) + "\n")
        n = len(lead_rows)
        for i in range(n):
            fh.write(",".join(list(lead_rows[i]) + [fmt_value(c, cols[c][i]) for c in columns]) + "\n")


def analyse_run(run_dir, args, quiet=False):
    """Analyse one run folder.  Returns (summary, output folder)."""
    run_dir = os.path.abspath(run_dir)
    if not os.path.isdir(run_dir):
        raise Stop("no such folder: %s" % run_dir)
    out_dir = os.path.abspath(args.out) if getattr(args, "out", None) else os.path.join(run_dir, "analysis_recheck")
    warnings = []
    manifest = None
    mpath = os.path.join(run_dir, "manifest.json")
    if os.path.isfile(mpath):
        try:
            with open(mpath) as fh:
                manifest = json.load(fh)
        except (ValueError, OSError) as err:
            warnings.append("manifest.json cannot be read (%s); it is not used" % err)
            manifest = None
    if not isinstance(manifest, dict):
        manifest = {} if manifest is None else {}
    has_manifest = os.path.isfile(mpath)

    # segment files
    segs = {"solute": [], "full": []}
    for f in sorted(os.listdir(run_dir)):
        m = SEG_RE.match(f)
        if m:
            segs[m.group(1)].append((int(m.group(2)), f))
    for k in segs:
        segs[k].sort()
    if not segs["solute"]:
        raise Stop("no solute_segNNNN.dcd in %s: nothing to analyse" % run_dir)

    who = identify_run(run_dir, args, manifest)
    warnings += who["warnings"]
    dt_s, dt_f = who["dt_s"], who["dt_f"]
    system = System(who["system_dir"])

    # topologies and guards on names and counts
    spdb = os.path.join(run_dir, "solute.pdb")
    if not os.path.isfile(spdb):
        raise Stop("solute.pdb is missing in %s" % run_dir)
    s_names, s_resn, _, _ = read_pdb(spdb)
    if len(s_names) != len(system.solute):
        raise Stop("solute.pdb has %d atoms, the list 'solute' of selections.json %d"
                   % (len(s_names), len(system.solute)))
    reader = SoluteReader(system)
    bad = system.name_mismatches(s_names, s_resn, lambda i: reader.pos.get(i))
    name_mismatches_solute = bad
    fpdb = os.path.join(run_dir, "full.pdb")
    have_full_pdb = os.path.isfile(fpdb)
    f_names = f_resn = None
    name_mismatches_full = []
    if have_full_pdb:
        f_names, f_resn, _, _ = read_pdb(fpdb)
        n_top = max([max(system.solute), max(system.water_O) if system.water_O else 0]) + 1
        if len(f_names) < n_top:
            raise Stop("full.pdb has %d atoms, but selections.json points to atom index %d"
                       % (len(f_names), n_top - 1))
        name_mismatches_full = system.name_mismatches(f_names, f_resn, lambda i: i)
    if name_mismatches_solute or name_mismatches_full:
        text = "atom names do not fit selections.json: " + "; ".join((name_mismatches_solute + name_mismatches_full)[:6])
        if getattr(args, "skip_name_check", False):
            warnings.append(text)
        else:
            raise Stop(text + "  (--skip-name-check turns this into a warning)")

    # ---------------- solute trajectory ----------------
    acc = reader.new_accumulator()
    sol_dcds = []
    overlap = {"suspected": False, "boundaries": []}
    prev_raw = None
    prev_name = None
    for idx, fname in segs["solute"]:
        d = Dcd(os.path.join(run_dir, fname))
        d.seg_index = idx
        sol_dcds.append(d)
        if not d.readable or d.n_frames == 0:
            continue
        if d.natoms != len(system.solute):
            raise Stop("%s has %d atoms per frame, solute.pdb %d" % (fname, d.natoms, len(system.solute)))
        raw = reader.read_file(d, acc, idx, keep_raw_ca=True)
        if prev_raw is not None:
            diff = np.abs(prev_raw - raw[0][None, :, :]).max(axis=(1, 2))
            k = int(diff.argmin())
            rep = bool(diff[k] < REPEAT_A)
            overlap["boundaries"].append({"segment_before": prev_name, "segment_after": fname,
                                          "repeated": rep,
                                          "frame_of_segment_before": k if rep else None,
                                          "smallest_largest_CA_difference_A": r6(float(diff[k]))})
            overlap["suspected"] = overlap["suspected"] or rep
        prev_raw, prev_name = raw, fname
    cols = SoluteReader.finish(acc)
    N = acc["n"]
    if N == 0:
        raise Stop("no complete production frame in %s (a run folder without a production frame is not analysed)"
                   % run_dir)
    seg_of_frame = list(acc["file_of_frame"])

    start_state = "C" if who["start"] in ("X", "N") else "R"
    opposite = "R" if start_state == "C" else "C"
    start_sign = 1 if who["start"] in ("X", "N") else -1
    state = [str(x) for x in cols["state"]]
    assigned = carry(state, ("C", "R"), start_state)
    cols["state"] = np.array(state)
    cols["state_assigned"] = np.array(assigned)
    sign = [int(x) for x in cols["n2_sign"]]
    sign_assigned = carry(sign, (1, -1), start_sign)
    cols["n2_sign_assigned"] = np.array(sign_assigned)

    # ---------------- settings of the water count of section 9D ----------------
    ring_place, w_settings = ring_place_settings(manifest)
    warnings += w_settings

    # ---------------- full trajectory ----------------
    full_rows = []
    full_aux = {"frames_without_box": 0}
    full_dcds = []
    if have_full_pdb:
        freader = FullReader(system, f_names, ring_place["site_atoms"], ring_place["water_cutoff_A"])
        for idx, fname in segs["full"]:
            d = Dcd(os.path.join(run_dir, fname))
            d.seg_index = idx
            full_dcds.append(d)
            if not d.readable or d.n_frames == 0:
                continue
            if d.natoms != len(f_names):
                raise Stop("%s has %d atoms per frame, full.pdb %d" % (fname, d.natoms, len(f_names)))
            freader.read_file(d, full_rows, full_aux, idx)
    else:
        for idx, fname in segs["full"]:
            d = Dcd(os.path.join(run_dir, fname))
            d.seg_index = idx
            full_dcds.append(d)
    NF = len(full_rows)
    fcols = {c: np.array([r[c] for r in full_rows]) for c in FULL_COLUMNS}
    if NF == 0:
        warnings.append("no full trajectory: no frame of full_segNNNN.dcd could be read%s; the quantities of "
                        "the full trajectory are missing and the full fraction of the fifth ligand does not exist"
                        % ("" if have_full_pdb else " (full.pdb is missing)"))

    # ---------------- guards of section 1 ----------------
    hc_s = header_check(sol_dcds, dt_s)
    hc_f = header_check(full_dcds, dt_f)
    for kind, hc in (("solute", hc_s), ("full", hc_f)):
        if hc["interval_matches"] is False:
            warnings.append("%s trajectory: the frame interval of a DCD header is not %g ps: %s"
                            % (kind, dt_s if kind == "solute" else dt_f, "; ".join(hc["details"][:3])))
        if hc["contiguous"] is False:
            warnings.append("%s trajectory: the DCD headers of the segments are not contiguous: %s"
                            % (kind, "; ".join(hc["details"][:3])))
    # manifest against the files
    mcheck = {"available": False, "consistent": None, "details": []}
    mseg = mget(manifest, "production", "segments")
    if isinstance(mseg, list):
        entries = [e for e in mseg if isinstance(e, dict) and isinstance(e.get("tag"), str)
                   and re.match(r"^seg\d+$", e["tag"])]
        if entries:
            mcheck["available"] = True
            mcheck["consistent"] = True
            read_s = {d.seg_index: d.n_frames for d in sol_dcds}
            read_f = {d.seg_index: d.n_frames for d in full_dcds}
            before = 0
            before_listed = True
            for e in sorted(entries, key=lambda e: int(e["tag"][3:])):
                i = int(e["tag"][3:])
                ns = mget(e, "frames", "solute")
                nf = mget(e, "frames", "full")
                got = read_s.get(i, 0)
                if is_number(ns) and int(ns) != got:
                    mcheck["consistent"] = False
                    mcheck["details"].append("%s: manifest %d solute frames, read %d" % (e["tag"], ns, got))
                if segs["full"] and is_number(nf) and int(nf) != read_f.get(i, 0):
                    mcheck["consistent"] = False
                    mcheck["details"].append("%s: manifest %d full frames, read %d" % (e["tag"], nf, read_f.get(i, 0)))
                frames_before = sum(n for j, n in read_s.items() if j < i)
                ts = e.get("time_start_ps")
                if is_number(ts) and abs(ts - frames_before * dt_s) > 1e-6 * max(1.0, dt_s):
                    mcheck["consistent"] = False
                    mcheck["details"].append("%s: manifest start %.6g ps, frames before it x interval = %.6g ps"
                                             % (e["tag"], ts, frames_before * dt_s))
            listed = {int(e["tag"][3:]) for e in entries}
            unlisted = [d.name for d in sol_dcds if d.seg_index not in listed and d.n_frames > 0]
            if unlisted:
                mcheck["details"].append("not in the list of the manifest (a segment in progress), not judged: "
                                         + ", ".join(unlisted))
            if mcheck["consistent"] is False:
                warnings.append("the segment list of the manifest does not agree with the files: "
                                + "; ".join(mcheck["details"][:4]))
    if overlap["suspected"]:
        b = [x for x in overlap["boundaries"] if x["repeated"]][0]
        warnings.append("repeated frame: the first frame of %s equals frame %d of %s (CA coordinates within %g A)"
                        % (b["segment_after"], b["frame_of_segment_before"], b["segment_before"], REPEAT_A))
    # system files against the manifest
    sys_match = {}
    for fn, path in (("selections.json", system.sel_path), ("reference.pdb", system.ref_path)):
        want = mget(manifest, "inputs", fn, "sha256")
        if isinstance(want, str):
            ok = sha256_file(path) == want
            sys_match[fn] = ok
            if not ok:
                warnings.append("%s of the system folder is not the file of the manifest (sha256 differs)" % fn)
        else:
            sys_match[fn] = None
    # whole protein
    largest_ca = acc["largest_neighbour_CA_distance_A"]
    if largest_ca > WHOLE_CA_A:
        warnings.append("the protein is not whole in the solute trajectory: neighbouring CA atoms up to %.3f A apart"
                        % largest_ca)
    # the two trajectories
    ratio = dt_f / dt_s
    cons = {"ratio_full_to_solute_interval": r6(ratio), "n_compared": 0, "max_abs_diff_d_Zn_N2_A": None}
    if NF and abs(ratio - round(ratio)) < 1e-9 and round(ratio) >= 1:
        k = int(round(ratio))
        diffs = [abs(fcols["d_Zn_N2_full_A"][j] - cols["d_Zn_N2_A"][k * (j + 1) - 1])
                 for j in range(NF) if k * (j + 1) - 1 < N]
        if diffs:
            cons["n_compared"] = len(diffs)
            cons["max_abs_diff_d_Zn_N2_A"] = r6(max(diffs))
            if max(diffs) > TWO_TRAJ_A:
                warnings.append("the full and the solute trajectory disagree: Zn-N2 differs by up to %.4f A "
                                "in frames of the same time" % max(diffs))

    # ---------------- time, completeness ----------------
    time_s = [(n + 1) * dt_s / 1000.0 for n in range(N)]
    time_f = [(j + 1) * dt_f / 1000.0 for j in range(NF)]
    seen_ns = N * dt_s / 1000.0
    planned = who["planned_ns"]
    n_planned = planned * 1000.0 / dt_s
    complete = bool(N >= n_planned * (1 - 1e-12))
    if who["role"] == "validation" and N > n_planned * (1 + 1e-12):
        warnings.append("validation run: %d solute frames, more than the %g planned; all frames read are used"
                        % (N, n_planned))

    # ---------------- statistics ----------------
    WINDOW_FRAMES = max(1, int(round(WINDOW_NS * 1000.0 / dt_s)))
    DWELL_FRAMES = max(1, int(round(DWELL_NS * 1000.0 / dt_s)))
    BLOCK_SOLUTE = max(1, int(round(BLOCK_NS * 1000.0 / dt_s)))
    BLOCK_FULL = max(1, int(round(BLOCK_NS * 1000.0 / dt_f)))
    stats = {}
    stat_cols = [c for c in SOLUTE_COLUMNS + (ND2_COLUMNS if system.has_nd2 else [])
                 if c.endswith("_A") or c.endswith("_deg") or c == "n2_volume"]
    for c in stat_cols:
        stats[c] = stat_block(cols[c], BLOCK_SOLUTE, "solute", c)
    for c in FULL_COLUMNS:
        if c.endswith("_A") or c in FULL_COUNTS:
            stats[c] = stat_block(fcols[c] if NF else [], BLOCK_FULL, "full", c)
    fractions = {}
    for f in FLAG_ORDER:
        if f in SOLUTE_FLAGS:
            fractions[f] = fraction_block(cols[f], BLOCK_SOLUTE)
        else:
            fractions[f] = fraction_block(fcols[f] if NF else [], BLOCK_FULL)
    if system.has_nd2:
        fractions["asn67_ND2_nearer_N1"] = fraction_block(cols["asn67_ND2_nearer_N1"], BLOCK_SOLUTE)
    lab_counts = {LABELLINGS[k]: int((cols["lig_rmsd_labelling"] == k).sum()) for k in range(4)}

    # ---------------- orientation (section 4) ----------------
    W = min(N, WINDOW_FRAMES)
    win = state[N - W:]
    cw = {s: win.count(s) for s in "CRI"}
    if 10 * cw[opposite] >= 9 * W:
        verdict = "CONVERTED"
    elif 10 * cw[start_state] >= 9 * W:
        verdict = "STAYED"
    else:
        verdict = "MIXED"
    visits = stretches(assigned)
    long_visits = []
    for st, first, n in visits:
        if n >= DWELL_FRAMES:
            strict = sum(1 for x in state[first:first + n] if x == st)
            long_visits.append({"state": st, "first_frame": first, "n_frames": n,
                                "t_first_ns": r6((first + 1) * dt_s / 1000.0),
                                "t_last_ns": r6((first + n) * dt_s / 1000.0), "long": True,
                                "fraction_strict": r6(strict / n)})
    seq = [start_state] + [v["state"] for v in long_visits]
    n_transitions = sum(1 for a, b in zip(seq[:-1], seq[1:]) if a != b)
    fp = next((v["first_frame"] for v in long_visits if v["state"] == opposite), None)
    h = N // 2
    orientation = {
        "start_state": start_state, "opposite_state": opposite, "midpoint_d1_A": MID_D1_A,
        "midpoint_d2_A": MID_D2_A, "window_frames": W, "window_ns": r6(W * dt_s / 1000.0),
        "window_complete": bool(N >= WINDOW_FRAMES), "count_C_window": cw["C"], "count_R_window": cw["R"],
        "count_I_window": cw["I"], "fraction_start_state_window": r6(cw[start_state] / W),
        "fraction_opposite_state_window": r6(cw[opposite] / W),
        "fraction_intermediate_window": r6(cw["I"] / W), "verdict": verdict, "provisional": not complete,
        "dwell_frames": DWELL_FRAMES,
        "first_passage_ns": r6((fp + 1) * dt_s / 1000.0) if fp is not None else None,
        "first_passage_frame": fp, "n_transitions": n_transitions, "n_visits": len(visits),
        "n_long_visits": len(long_visits), "long_visits": long_visits, "state_first_frame": state[0],
        "count_C": state.count("C"), "count_R": state.count("R"), "count_I": state.count("I"),
        "fraction_C": r6(state.count("C") / N), "fraction_R": r6(state.count("R") / N),
        "fraction_I": r6(state.count("I") / N), "n_last_half": h,
        "count_C_last_half": state[N - h:].count("C") if h else 0,
        "count_R_last_half": state[N - h:].count("R") if h else 0,
        "count_I_last_half": state[N - h:].count("I") if h else 0}

    # ---------------- Asn67 ND2 ----------------
    if system.has_nd2:
        fl = [int(x) for x in cols["asn67_ND2_nearer_N1"]]
        by = {}
        for s_ in "CRI":
            idx = [i for i in range(N) if state[i] == s_]
            by[s_] = {"n": len(idx), "count_ND2_nearer": sum(fl[i] for i in idx)}
        asn = {"available": True, "note": "reported only; the states use OD1", "n": N,
               "count_ND2_nearer": sum(fl), "fraction_ND2_nearer": r6(sum(fl) / N), "by_state": by}
    else:
        asn = {"available": False, "note": "selections.json has no key asn67_ND2"}

    # ---------------- configuration at N2 (section 9A) ----------------
    def sign_counts(sub):
        return sub.count(1), sub.count(-1), sub.count(0)
    cp, cn, cz = sign_counts(sign)
    wp, wn, wz = sign_counts(sign[N - W:])
    by_state = {}
    for s_ in "CRI":
        sub = [sign[i] for i in range(N) if state[i] == s_]
        a_, b_, c_ = sign_counts(sub)
        by_state[s_] = {"n": len(sub), "count_positive": a_, "count_negative": b_, "count_zero": c_,
                        "fraction_positive": r6(a_ / len(sub)) if sub else None,
                        "fraction_negative": r6(b_ / len(sub)) if sub else None}
    sstr = stretches(sign_assigned)
    long_str = [{"sign": int(v), "first_frame": f, "n_frames": n, "t_first_ns": r6((f + 1) * dt_s / 1000.0),
                 "t_last_ns": r6((f + n) * dt_s / 1000.0), "long": True}
                for v, f, n in sstr if n >= DWELL_FRAMES]
    changes = []
    prev = start_sign
    for ls in long_str:
        if ls["sign"] != prev:
            changes.append({"frame": ls["first_frame"], "time_ns": ls["t_first_ns"], "from_sign": prev,
                            "to_sign": ls["sign"]})
        prev = ls["sign"]
    n2conf = {
        "order": "Zn, HN21, S", "sign_of_start": start_sign, "n": N, "count_positive": cp,
        "count_negative": cn, "count_zero": cz, "fraction_positive": r6(cp / N),
        "fraction_negative": r6(cn / N), "window_frames": W, "count_positive_window": wp,
        "count_negative_window": wn, "count_zero_window": wz, "fraction_positive_window": r6(wp / W),
        "fraction_negative_window": r6(wn / W), "by_state": by_state, "sign_first_frame": sign[0],
        "sign_last_frame": sign[-1], "dwell_frames": DWELL_FRAMES, "n_stretches": len(sstr),
        "n_long_stretches": len(long_str), "long_stretches": long_str,
        "n_lasting_sign_changes": len(changes), "lasting_sign_changes": changes,
        "first_lasting_sign_change_ns": changes[0]["time_ns"] if changes else None}

    # ---------------- hydrogen bond before and after the first passage (9B) ----------------
    hb = [int(x) for x in cols["hb_N2H_T199OG1"]]
    if fp is not None:
        cb = sum(hb[:fp])
        ca_ = sum(hb[fp:])
        hbfp = {"available": True, "flag": "hb_N2H_T199OG1", "verdict": verdict, "first_passage_frame": fp,
                "first_passage_ns": r6((fp + 1) * dt_s / 1000.0), "n_before": fp, "count_before": cb,
                "fraction_before": r6(cb / fp) if fp else None, "n_after": N - fp, "count_after": ca_,
                "fraction_after": r6(ca_ / (N - fp))}
    else:
        hbfp = {"available": False, "flag": "hb_N2H_T199OG1", "verdict": verdict, "first_passage_frame": None,
                "first_passage_ns": None, "n_before": None, "count_before": None, "fraction_before": None,
                "n_after": None, "count_after": None, "fraction_after": None}

    # ---------------- first production frame against the start (9E) ----------------
    sop = {"note": "solute frame 0 of the production against the start as built",
           "time_ns": r6(dt_s / 1000.0), "state_of_start": start_state, "state_first_frame": state[0],
           "state_differs": state[0] != start_state, "state_first_frame_is_opposite": state[0] == opposite,
           "n2_sign_of_start": start_sign, "n2_sign_first_frame": sign[0],
           "n2_sign_differs": sign[0] != start_sign, "n2_volume_first_frame": r6(cols["n2_volume"][0]),
           "lig_rmsd_first_frame_A": r6(cols["lig_rmsd_A"][0])}

    # ---------------- RMSF (section 9) ----------------
    mean_dev = acc["sum_dev"] / N
    rmsf = np.sqrt(np.clip(acc["sum_dev2"] / N - (mean_dev ** 2).sum(axis=1), 0.0, None))
    kmax = int(rmsf.argmax())
    rmsf_block = {"file": "rmsf_ca.csv", "n_residues": len(rmsf), "mean_A": r6(rmsf.mean()),
                  "max_A": r6(rmsf[kmax]), "max_residue_number_deposited": system.resnum_ca[kmax]}

    # ---------------- validation numbers and gate G-V (sections 5 and 11) ----------------
    st_l = stats["lig_rmsd_A"]
    rule1 = bool(st_l["mean_last_half"] is not None and st_l["mean_last_half"] > RULE1_A)
    validation = {"lig_rmsd_mean_A": st_l["mean"], "lig_rmsd_sd_A": st_l["sd"],
                  "lig_rmsd_mean_last_half_A": st_l["mean_last_half"], "rule1_threshold_A": RULE1_A,
                  "rule1_left_crystal": rule1,
                  "rule1_applies_to_this_run": bool(who["start"] == "X" and who["role"] == "production"),
                  "hb_N2H_T199OG1_fraction": fractions["hb_N2H_T199OG1"]["fraction"]}
    is_val = who["role"] == "validation"
    gv1 = {"mean_A": st_l["mean"], "limit_max_A": GV1_MAX_A, "passed": bool(st_l["mean"] <= GV1_MAX_A)}
    fr = fractions["hb_N2H_T199OG1"]
    gv2 = {"fraction": fr["fraction"], "count": fr["count"], "n": fr["n"], "limit_min": 0.8,
           "passed": bool(5 * fr["count"] >= 4 * fr["n"])}
    crit = {}
    for col, win_ in [("d_Zn_N2_A", GV3A_WINDOW_A), ("d_Zn_H94NE2_A", GV3B_WINDOW_A),
                      ("d_Zn_H96NE2_A", GV3B_WINDOW_A), ("d_Zn_H119ND1_A", GV3B_WINDOW_A)]:
        m_ = stats[col]["mean"]
        crit[col] = {"mean_A": m_, "window_A": [win_[0], win_[1]], "passed": bool(win_[0] <= m_ <= win_[1])}
    for col, ref in GV3C_REF_DEG:
        m_ = stats[col]["mean"]
        lo, hi = round(ref - GV3C_TOL_DEG, 6), round(ref + GV3C_TOL_DEG, 6)
        crit[col] = {"mean_deg": m_, "reference_6RL9_deg": ref, "tolerance_deg": GV3C_TOL_DEG,
                     "deviation_deg": r6(m_ - ref), "window_deg": [lo, hi], "passed": bool(lo <= m_ <= hi)}
    fs, fa = fractions["fifth_sulfonylO"], fractions["fifth_any"]
    sol_ok = bool(20 * fs["count"] <= fs["n"])
    full_ok = bool(NF >= 1 and 20 * fa["count"] <= NF)
    crit["fifth_ligand"] = {
        "limit_max_fraction": 0.05,
        "solute_trajectory": {"flag": "fifth_sulfonylO", "count": fs["count"], "n": fs["n"],
                              "fraction": fs["fraction"], "passed": sol_ok},
        "full_trajectory": {"flag": "fifth_any", "available": bool(NF >= 1),
                            "count": fa["count"] if NF else None, "n": NF,
                            "fraction": fa["fraction"] if NF else None, "passed": full_ok},
        "passed": bool(sol_ok and full_ok)}
    gv3 = {"criteria": crit, "passed": all(v["passed"] for v in crit.values())}
    c4 = {}
    for col, ref, tol in (("ang_C4_S_N2_deg", GV4A_REF_DEG, GV4A_TOL_DEG), ("tau_fold_deg", GV4B_REF_DEG, GV4B_TOL_DEG)):
        m_ = stats[col]["mean"]
        lo, hi = round(ref - tol, 6), round(ref + tol, 6)
        c4[col] = {"mean_deg": m_, "sd_deg": stats[col]["sd"], "reference_deg": ref, "tolerance_deg": tol,
                   "deviation_deg": r6(m_ - ref), "window_deg": [lo, hi], "passed": bool(lo <= m_ <= hi)}
    gv4 = {"criteria": c4, "passed": all(v["passed"] for v in c4.values())}
    gv_passed = bool(gv1["passed"] and gv2["passed"] and gv3["passed"] and gv4["passed"])
    glu = {"note": "reported beside the gate, no criterion (sections 9C and 11)", "distance_max_A": FIFTH_A}
    for col, flag in (("d_Zn_E106OE1_A", "fifth_E106OE1"), ("d_Zn_E106OE2_A", "fifth_E106OE2")):
        glu[col] = {"mean_A": stats[col]["mean"], "min_A": stats[col]["min"],
                    "count_within": fractions[flag]["count"], "n": fractions[flag]["n"],
                    "fraction_within": fractions[flag]["fraction"]}
    gate = {"applies_to_this_run": is_val,
            "note": "gate G-V of md/PLAN.md; a verdict only for the validation run X0",
            "GV1_lig_rmsd": gv1, "GV2_hbond": gv2, "GV3_zinc": gv3, "GV4_ligand_internal": gv4,
            "passed": gv_passed, "verdict": ("PASSED" if gv_passed else "FAILED") if is_val else None,
            "provisional": not complete, "reported_zn_glu106": glu}

    # ---------------- equilibration (section 9E) ----------------
    os.makedirs(out_dir, exist_ok=True)
    equil = analyse_equilibration(run_dir, out_dir, manifest, system, reader, f_names if have_full_pdb else None,
                                  dt_s, dt_f, start_state, opposite, start_sign, warnings, ring_place)

    # ---------------- write the series ----------------
    columns = SOLUTE_COLUMNS + (ND2_COLUMNS if system.has_nd2 else [])
    lead = [("%d" % n, "%d" % seg_of_frame[n], "%.3f" % time_s[n]) for n in range(N)]
    write_solute_csv(os.path.join(out_dir, "timeseries_solute.csv"), cols, ["frame", "segment", "time_ns"],
                     lead, columns)
    lead = [("%d" % j, "%d" % full_rows[j]["_tag"], "%.3f" % time_f[j]) for j in range(NF)]
    write_solute_csv(os.path.join(out_dir, "timeseries_full.csv"), fcols, ["frame", "segment", "time_ns"],
                     lead, FULL_COLUMNS)
    with open(os.path.join(out_dir, "rmsf_ca.csv"), "w") as fh:
        fh.write("residue_number_deposited,residue_name,rmsf_A\n")
        for i in range(len(rmsf)):
            fh.write("%d,%s,%.4f\n" % (system.resnum_ca[i], system.ref_ca_resn[i], rmsf[i]))

    # frames near a threshold, for the comparison of two programs (section 14)
    aux = {"what": "frames in which this program finds a value within the tolerance of section 14 of a "
                   "threshold; used by --compare to tell a threshold frame from a disagreement",
           "labelling_tie_frames": [int(i) for i in np.where(cols["_labelling_gap_A"] < 0.001)[0]],
           "full_frames_water_N1_at_3p5A": [j for j in range(NF) if full_rows[j]["_edge_n1"]],
           "full_frames_water_ring_wide_at_3p5A": [j for j in range(NF) if full_rows[j]["_edge_wide"]],
           "full_frames_water_ring_place_at_cutoff": [j for j in range(NF) if full_rows[j]["_edge_place"]],
           "water_ring_place_cutoff_A": ring_place["water_cutoff_A"]}
    aux.update(equil.pop("_aux", {}))
    write_json(os.path.join(out_dir, "recheck_aux.json"), aux)

    # ---------------- checksums ----------------
    inputs = {"selections.json": sha256_file(system.sel_path), "reference.pdb": sha256_file(system.ref_path),
              "solute.pdb": sha256_file(spdb), "full.pdb": sha256_file(fpdb) if have_full_pdb else None,
              "manifest.json": sha256_file(mpath) if has_manifest else None}
    dcd_names = sorted(f for f in os.listdir(run_dir) if SEG_RE.match(f) or EQ_FILE_RE.match(f))
    for f in dcd_names:
        if getattr(args, "no_hash_trajectories", False):
            inputs[f] = "not computed (--no-hash-trajectories); %d bytes" % os.path.getsize(os.path.join(run_dir, f))
        else:
            inputs[f] = sha256_file(os.path.join(run_dir, f))

    # ---------------- the line of the campaign table ----------------
    def mean_of(col, key="mean"):
        return stats[col][key] if col in stats else None
    ang_means = [stats[c]["mean"] for c, _ in GV3C_REF_DEG]
    ang_dev = [abs(stats[c]["mean"] - ref) for c, ref in GV3C_REF_DEG]
    row = {
        "run": who["label"], "model": who["model"], "start": who["start"], "replica": who["replica"],
        "seen_ns": r6(seen_ns), "planned_ns": r6(planned), "complete": complete, "verdict": verdict,
        "provisional": not complete,
        "frac_start_state_last20ns": orientation["fraction_start_state_window"],
        "frac_opposite_state_last20ns": orientation["fraction_opposite_state_window"],
        "frac_intermediate_last20ns": orientation["fraction_intermediate_window"],
        "first_passage_ns": orientation["first_passage_ns"], "n_transitions": n_transitions,
        "d1_N1_F131CZ_mean_last_half_A": mean_of("d1_N1_F131CZ_A", "mean_last_half"),
        "d2_N1_N67OD1_mean_last_half_A": mean_of("d2_N1_N67OD1_A", "mean_last_half"),
        "d_N1_N67ND2_mean_last_half_A": mean_of("d_N1_N67ND2_A", "mean_last_half") if system.has_nd2 else None,
        "asn67_ND2_nearer_N1_fraction": fractions["asn67_ND2_nearer_N1"]["fraction"] if system.has_nd2 else None,
        "lig_rmsd_mean_A": st_l["mean"], "lig_rmsd_sd_A": st_l["sd"],
        "lig_rmsd_mean_last_half_A": st_l["mean_last_half"], "lig_rmsd_sd_last_half_A": st_l["sd_last_half"],
        "rule1_left_crystal": rule1,
        "hb_N2H_T199OG1_fraction": fractions["hb_N2H_T199OG1"]["fraction"],
        "hb_T199NH_Osulfonyl_fraction": fractions["hb_T199NH_Osulfonyl"]["fraction"],
        "hb_T199OG1H_E106OE_fraction": fractions["hb_T199OG1H_E106OE"]["fraction"],
        "d_Zn_N2_mean_A": mean_of("d_Zn_N2_A"), "d_Zn_N2_sd_A": mean_of("d_Zn_N2_A", "sd"),
        "d_Zn_H94NE2_mean_A": mean_of("d_Zn_H94NE2_A"), "d_Zn_H96NE2_mean_A": mean_of("d_Zn_H96NE2_A"),
        "d_Zn_H119ND1_mean_A": mean_of("d_Zn_H119ND1_A"),
        "ang_N_Zn_N_smallest_mean_deg": r6(min(ang_means)), "ang_N_Zn_N_largest_mean_deg": r6(max(ang_means)),
        "ang_N_Zn_N_largest_abs_deviation_deg": r6(max(ang_dev)),
        "d_Zn_Osulfonyl_min_mean_A": mean_of("d_Zn_Osulfonyl_min_A"),
        "fifth_sulfonylO_fraction": fractions["fifth_sulfonylO"]["fraction"],
        "d_Zn_Owater_min_mean_A": mean_of("d_Zn_Owater_min_A"),
        "d_Zn_Owater_min_min_A": mean_of("d_Zn_Owater_min_A", "min"),
        "fifth_any_fraction": fractions["fifth_any"]["fraction"],
        "n_water_N1_3p5A_mean": mean_of("n_water_N1_3p5A"),
        "tau_fold_mean_deg": mean_of("tau_fold_deg"), "tau_fold_sd_deg": mean_of("tau_fold_deg", "sd"),
        "ang_C4_S_N2_mean_deg": mean_of("ang_C4_S_N2_deg"), "ang_C4_S_N2_sd_deg": mean_of("ang_C4_S_N2_deg", "sd"),
        "d_S_N2_mean_A": mean_of("d_S_N2_A"), "d_S_N2_sd_A": mean_of("d_S_N2_A", "sd"),
        "ca_rmsd_mean_A": mean_of("ca_rmsd_A"), "ca_rmsd_mean_last_half_A": mean_of("ca_rmsd_A", "mean_last_half"),
        "rmsf_ca_max_A": rmsf_block["max_A"], "rmsf_ca_max_residue": rmsf_block["max_residue_number_deposited"],
        "gate_GV_passed": gv_passed, "n_warnings": None,
        "n2_fraction_positive": n2conf["fraction_positive"],
        "n2_fraction_positive_last20ns": n2conf["fraction_positive_window"],
        "n2_n_lasting_sign_changes": n2conf["n_lasting_sign_changes"],
        "n2_first_lasting_sign_change_ns": n2conf["first_lasting_sign_change_ns"],
        "hb_N2H_T199OG1_fraction_before_first_passage": hbfp["fraction_before"],
        "hb_N2H_T199OG1_fraction_after_first_passage": hbfp["fraction_after"],
        "d_Zn_E106OE1_mean_A": mean_of("d_Zn_E106OE1_A"), "d_Zn_E106OE1_min_A": mean_of("d_Zn_E106OE1_A", "min"),
        "fifth_E106OE1_fraction": fractions["fifth_E106OE1"]["fraction"],
        "d_Zn_E106OE2_mean_A": mean_of("d_Zn_E106OE2_A"), "d_Zn_E106OE2_min_A": mean_of("d_Zn_E106OE2_A", "min"),
        "fifth_E106OE2_fraction": fractions["fifth_E106OE2"]["fraction"],
        "n_water_ring_wide_3p5A_mean": mean_of("n_water_ring_wide_3p5A"),
        "d_water_ring_wide_min_min_A": mean_of("d_water_ring_wide_min_A", "min"),
        "state_first_frame": state[0], "state_first_frame_differs_from_start": sop["state_differs"],
        "n2_sign_first_frame": sign[0],
        "equilibration_frames_solute": equil.get("n_frames_solute", 0) if equil.get("available") else 0,
        "n_water_ring_place_mean": mean_of("n_water_ring_place"),
        "d_water_ring_place_min_min_A": mean_of("d_water_ring_place_min_A", "min")}
    row["n_warnings"] = len(warnings)
    row = {k: row[k] for k in TABLE_ROW_KEYS}

    summary = {
        "spec": "md/analysis/SPEC.md version " + SPEC_VERSION,
        "script": {"name": PROGRAM, "sha256": sha256_file(os.path.abspath(__file__)),
                   "what": "second, independent implementation of the specification"},
        "software": {"python": platform.python_version(), "numpy": np.__version__,
                     "MDAnalysis": "not used", "matplotlib": "not used"},
        "random_numbers": "none",
        "units": "distances in Angstrom, angles in degrees, times in ns, frame intervals in ps",
        "run": {
            "run_dir": run_dir, "label": who["label"], "model": who["model"], "start": who["start"],
            "replica": who["replica"], "role": who["role"], "system_dir": who["system_dir"],
            "planned_ns": r6(planned), "seen_ns": r6(seen_ns), "complete": complete,
            "solute_interval_ps": dt_s, "full_interval_ps": dt_f, "n_frames_solute": N,
            "n_frames_full": NF, "seen_full_ns": r6(NF * dt_f / 1000.0),
            "n_atoms_solute": len(system.solute), "n_CA": len(system.ca),
            "segments_solute": [d.info(d.seg_index) for d in sol_dcds],
            "segments_full": [d.info(d.seg_index) for d in full_dcds],
            "header_check_solute": hc_s, "header_check_full": hc_f, "segment_overlap_check": overlap,
            "manifest_segment_check": mcheck, "sources": who["sources"],
            "system_files_match_manifest": sys_match,
            "imaging_solute": {k: acc[k] for k in ("frames_with_box", "frames_without_box",
                                                   "frames_zinc_moved", "frames_ligand_moved")},
            "full_trajectory": {"n_water_O": len(system.water_O),
                                "n_protein_O": int(len(freader.prot_O)) if have_full_pdb else None,
                                "frames_without_box": full_aux["frames_without_box"],
                                "name_mismatches": name_mismatches_full},
            "largest_neighbour_CA_distance_A": r6(largest_ca), "consistency_full_vs_solute": cons},
        "inputs_sha256": inputs,
        "constants": {
            "midpoint_d1_A": MID_D1_A, "midpoint_d2_A": MID_D2_A, "hbond_distance_max_A": HB_DIST_A,
            "hbond_angle_min_deg": HB_ANGLE_DEG, "fifth_ligand_distance_max_A": FIFTH_A,
            "water_distance_max_A": WATER_A, "ring_place_default_site_atoms": RING_PLACE_ATOMS,
            "ring_place_default_cutoff_A": RING_PLACE_CUTOFF_A, "ring_wide_site_atoms": RING_WIDE_ATOMS,
            "ring_wide_cutoff_A": WATER_A, "driver_nearest_tolerance_A": DRIVER_NEAREST_TOL_A,
            "window_ns": WINDOW_NS, "dwell_ns": DWELL_NS, "block_ns": BLOCK_NS,
            "window_frames": WINDOW_FRAMES, "dwell_frames": DWELL_FRAMES,
            "block_frames_solute": BLOCK_SOLUTE, "block_frames_full": BLOCK_FULL, "rule1_threshold_A": RULE1_A,
            "validation_planned_ns": VALIDATION_NS, "akma_time_unit_ps": AKMA_PS,
            "equilibration_stage_table_ps": EQ_CONSTANT_TABLE},
        "stats": stats, "fractions": fractions, "lig_rmsd_labelling_counts": lab_counts,
        "orientation": orientation, "asn67_ND2": asn, "n2_configuration": n2conf,
        "hbond_first_passage": hbfp, "water_ring_place": ring_place, "start_of_production": sop,
        "equilibration": equil,
        "validation": validation, "gate_GV": gate, "rmsf": rmsf_block, "warnings": warnings,
        "table_row": row}
    write_json(os.path.join(out_dir, "summary.json"), summary)
    if not quiet:
        print("%s: %d solute and %d full frames, %.3f of %.3f ns%s; verdict %s%s; ligand RMSD %.3f A; "
              "gate G-V numbers %s; %d warning(s); written to %s"
              % (who["label"], N, NF, seen_ns, planned, "" if complete else " (not complete)", verdict,
                 " (provisional)" if not complete else "", st_l["mean"],
                 ("PASSED" if gv_passed else "FAILED") if is_val else ("met" if gv_passed else "not met") + ", no verdict",
                 len(warnings), out_dir))
        for w_ in warnings:
            print("   warning: " + w_)
    return summary, out_dir


# ----------------------------------------------------------------------------------------------
# equilibration frames (section 9E)
# ----------------------------------------------------------------------------------------------
def analyse_equilibration(run_dir, out_dir, manifest, system, reader, full_names, dt_s, dt_f,
                          start_state, opposite, start_sign, warnings, ring_place):
    files = {"solute": {}, "full": {}}
    for f in sorted(os.listdir(run_dir)):
        m = EQ_FILE_RE.match(f)
        if m:
            files[m.group(1)][m.group(2)] = f
    listed_s = sorted(files["solute"].values())
    listed_f = sorted(files["full"].values())
    csv_s = os.path.join(out_dir, "timeseries_equil_solute.csv")
    csv_f = os.path.join(out_dir, "timeseries_equil_full.csv")

    def remove_old(paths):
        for p in paths:
            if os.path.isfile(p):
                os.remove(p)
    if not listed_s and not listed_f:
        remove_old([csv_s, csv_f])
        return {"available": False, "note": "no equilibration files in the run folder",
                "files_solute": [], "files_full": []}
    # step table
    stages = mget(manifest, "equilibration", "stages")
    steps = []
    if isinstance(stages, list):
        source = "manifest equilibration.stages"
        for e in stages:
            if not isinstance(e, dict) or not isinstance(e.get("tag"), str):
                continue
            m = EQ_TAG_RE.match(e["tag"])
            if m and is_number(e.get("time_start_ps")) and is_number(e.get("time_end_ps")):
                steps.append({"tag": e["tag"], "stage": int(m.group(1)), "extension": int(m.group(2) or 0),
                              "S": float(e["time_start_ps"]), "E": float(e["time_end_ps"]), "entry": e})
    else:
        source = "protocol (constants of the specification)"
        for nn, s_, e_ in EQ_CONSTANT_TABLE:
            steps.append({"tag": "eq%02d" % nn, "stage": nn, "extension": 0, "S": s_, "E": e_, "entry": None})
    steps = list({st["tag"]: st for st in steps}.values())        # a tag listed twice: the last entry counts
    steps.sort(key=lambda x: (x["stage"], x["extension"]))
    known = {st["tag"] for st in steps}
    if not steps:
        warnings.append("equilibration: the step table is empty; no equilibration file is read")
    unknown = sorted((set(files["solute"]) | set(files["full"])) - known)
    for tag in unknown:
        if steps:
            warnings.append("equilibration: the files of step %s (%s) belong to no step of the step table "
                            "(%s) and are not read" % (tag, ", ".join(files[k][tag] for k in ("solute", "full")
                                                                      if tag in files[k]), source))
    T_eq = max([st["E"] for st in steps]) if steps else 0.0
    have_s = [st for st in steps if st["tag"] in files["solute"]]
    have_f = [st for st in steps if st["tag"] in files["full"]]
    if not have_s and not have_f:
        remove_old([csv_s, csv_f])
        return {"available": False, "note": "equilibration files are in the run folder, but none belongs to a "
                                            "step of the step table", "files_solute": listed_s,
                "files_full": listed_f, "stage_table_source": source}
    if have_f and full_names is None:
        warnings.append("equilibration: full_eqNN.dcd without full.pdb; the full frames are not read")

    def open_file(kind, st, dt):
        d = Dcd(os.path.join(run_dir, files[kind][st["tag"]]))
        expected = int(math.floor((st["E"] - st["S"]) * (1 + 1e-9) / dt))
        if not d.readable:
            warnings.append("%s: cannot be read as a DCD file (%s); taken as holding no frame" % (d.name, d.error))
        want_atoms = len(system.solute) if kind == "solute" else (len(full_names) if full_names else None)
        if d.readable and want_atoms is not None and d.natoms != want_atoms:
            raise Stop("%s has %d atoms per frame, %d expected" % (d.name, d.natoms, want_atoms))
        if d.n_frames != expected:
            warnings.append("%s: %d frames read, %d expected for step %s (%g to %g ps)"
                            % (d.name, d.n_frames, expected, st["tag"], st["S"], st["E"]))
        if d.readable and d.n_frames > 0:
            if abs(d.interval_ps - dt) > 0.001 * dt:
                warnings.append("%s: the frame interval of the DCD header is %.6g ps, not %g ps"
                                % (d.name, d.interval_ps, dt))
            th = d.first_frame_ps
            if abs(th - (st["S"] + dt)) > 0.001 * dt and abs(th - dt) > 0.001 * dt:
                warnings.append("%s: the DCD header puts the first frame at %.6g ps; expected %.6g ps (or %g ps)"
                                % (d.name, th, st["S"] + dt, dt))
        info = d.info()
        info.update({"tag": st["tag"], "stage": st["stage"], "extension": st["extension"],
                     "n_frames_expected": expected})
        if d.readable and d.n_frames > 0:
            info["header_first_frame_ps"] = r6(d.first_frame_ps)
        return d, info

    # solute frames
    acc = reader.new_accumulator()
    per_step = {}
    lead = []
    frame = 0
    for st in have_s:
        d, info = open_file("solute", st, dt_s)
        first = frame
        if d.readable and d.n_frames:
            reader.read_file(d, acc, st["tag"])
            for k in range(d.n_frames):
                t = (st["S"] + (k + 1) * dt_s - T_eq) / 1000.0
                lead.append(("%d" % frame, "%d" % st["stage"], "%d" % st["extension"], "%.6f" % t, t))
                frame += 1
        per_step[st["tag"]] = {"solute": info, "first": first, "n": frame - first}
    cols = SoluteReader.finish(acc)
    NS = frame
    columns = [c for c in SOLUTE_COLUMNS + (ND2_COLUMNS if system.has_nd2 else [])
               if c not in ("state_assigned", "n2_sign_assigned")]
    written = []
    if have_s:
        if NS:
            cols["state"] = np.array([str(x) for x in cols["state"]])
        else:
            cols = {c: [] for c in columns}
        write_solute_csv(csv_s, cols, ["frame", "stage", "extension", "time_ns"], [x[:4] for x in lead], columns)
        written.append("timeseries_equil_solute.csv")
    else:
        remove_old([csv_s])

    # full frames
    frows = []
    flead = []
    per_step_full = {}
    eaux = {"equil_full_frames_water_N1_at_3p5A": [], "equil_full_frames_water_ring_wide_at_3p5A": [],
            "equil_full_frames_water_ring_place_at_cutoff": [], "equil_driver_check": {}}
    faux = {"frames_without_box": 0}
    if have_f and full_names is not None:
        freader = FullReader(system, full_names, ring_place["site_atoms"], ring_place["water_cutoff_A"])
        frame = 0
        for st in have_f:
            d, info = open_file("full", st, dt_f)
            first = frame
            if d.readable and d.n_frames:
                freader.read_file(d, frows, faux, st["tag"])
                for k in range(d.n_frames):
                    t = (st["S"] + (k + 1) * dt_f - T_eq) / 1000.0
                    flead.append(("%d" % frame, "%d" % st["stage"], "%d" % st["extension"], "%.6f" % t, t))
                    frame += 1
            per_step_full[st["tag"]] = {"full": info, "first": first, "n": frame - first}
        fc = {c: [r[c] for r in frows] for c in FULL_COLUMNS}
        write_solute_csv(csv_f, fc, ["frame", "stage", "extension", "time_ns"], [x[:4] for x in flead], FULL_COLUMNS)
        written.append("timeseries_equil_full.csv")
        eaux["equil_full_frames_water_N1_at_3p5A"] = [j for j, r in enumerate(frows) if r["_edge_n1"]]
        eaux["equil_full_frames_water_ring_wide_at_3p5A"] = [j for j, r in enumerate(frows) if r["_edge_wide"]]
        eaux["equil_full_frames_water_ring_place_at_cutoff"] = [j for j, r in enumerate(frows) if r["_edge_place"]]
    else:
        remove_old([csv_f])
    NFe = len(frows)

    def frame_entry(i):
        return {"frame": i, "stage": int(lead[i][1]), "extension": int(lead[i][2]), "time_ns": r6(lead[i][4]),
                "state": str(cols["state"][i]), "n2_sign": int(cols["n2_sign"][i])}
    stage_list = []
    for st in steps:
        tag = st["tag"]
        if tag not in files["solute"] and tag not in files["full"]:
            continue
        ps = per_step.get(tag)
        pf = per_step_full.get(tag)
        e = {"stage": st["stage"], "extension": st["extension"], "tag": tag, "time_start_ps": r6(st["S"]),
             "time_end_ps": r6(st["E"]), "time_start_ns": r6((st["S"] - T_eq) / 1000.0),
             "time_end_ns": r6((st["E"] - T_eq) / 1000.0), "solute": ps["solute"] if ps else None,
             "full": pf["full"] if pf else None}
        if ps and ps["n"]:
            sl = slice(ps["first"], ps["first"] + ps["n"])
            stt = list(cols["state"][sl])
            sg = [int(x) for x in cols["n2_sign"][sl]]
            last = ps["first"] + ps["n"] - 1
            e.update({"count_C": stt.count("C"), "count_R": stt.count("R"), "count_I": stt.count("I"),
                      "count_n2_positive": sg.count(1), "count_n2_negative": sg.count(-1),
                      "last_solute_frame": {"frame": last, "time_ns": r6(lead[last][4]),
                                            "state": str(cols["state"][last]),
                                            "lig_rmsd_A": r6(cols["lig_rmsd_A"][last]),
                                            "n2_volume": r6(cols["n2_volume"][last]),
                                            "n2_sign": int(cols["n2_sign"][last])}})
        else:
            e.update({"count_C": 0, "count_R": 0, "count_I": 0, "count_n2_positive": 0, "count_n2_negative": 0,
                      "last_solute_frame": None})
        last_full = None
        m_full = pf["n"] if pf else 0           # frames read from the full file of the step
        if m_full:
            lf = pf["first"] + m_full - 1
            last_full = {"frame": lf, "time_ns": r6(flead[lf][4]),
                         "n_water_ring_place": int(frows[lf]["n_water_ring_place"]),
                         "d_water_ring_place_min_A": r6(frows[lf]["d_water_ring_place_min_A"]),
                         "n_water_ring_wide_3p5A": int(frows[lf]["n_water_ring_wide_3p5A"]),
                         "d_water_ring_wide_min_A": r6(frows[lf]["d_water_ring_wide_min_A"])}
        e["last_full_frame"] = last_full
        # the driver's own count: written for a step whose entry of the manifest holds
        # hydration_check.n_water_O as a JSON integer; compared with the first count of section 9D
        # in the last full frame of the step if that frame is written at the end of the step
        hc = mget(st["entry"], "hydration_check") if st["entry"] else None
        dh = None
        if isinstance(hc, dict) and is_json_integer(hc.get("n_water_O")):
            near = mget(hc, "nearest_water_O", "distance_A")
            near_is_number = is_number(near)
            dh = {"n_water_O": int(hc["n_water_O"]), "nearest_water_O_A": r6(near) if near_is_number else near,
                  "applies": hc.get("applies"), "decision": hc.get("decision"), "compared": False,
                  "count_agrees": None, "nearest_agrees": None}
            if m_full >= 1 and abs(st["S"] + m_full * dt_f - st["E"]) <= STEP_END_TOL_PS:
                d_here = float(frows[lf]["d_water_ring_place_min_A"])
                dh["compared"] = True
                dh["count_agrees"] = bool(int(frows[lf]["n_water_ring_place"]) == dh["n_water_O"])
                if near_is_number and math.isfinite(d_here):
                    dh["nearest_agrees"] = bool(abs(d_here - near) <= DRIVER_NEAREST_TOL_A)
                eaux["equil_driver_check"][tag] = {
                    "full_frame": lf, "water_oxygen_at_cutoff": bool(frows[lf]["_edge_place"]),
                    "abs_difference_from_driver_A": abs(d_here - near) if near_is_number and math.isfinite(d_here) else None}
                if dh["count_agrees"] is False or dh["nearest_agrees"] is False:
                    warnings.append("equilibration: the driver's count of water in the place of the ring at the "
                                    "end of %s (%d, nearest %s A) does not agree with the full frame of that "
                                    "time (n_water_ring_place %d, nearest %.4f A)"
                                    % (tag, dh["n_water_O"], near, int(frows[lf]["n_water_ring_place"]), d_here))
        e["driver_hydration_check"] = dh
        stage_list.append(e)

    out = {"available": True, "note": "equilibration frames: a series of their own, times zero or negative; "
                                      "they enter no statistic of the production",
           "files_solute": listed_s, "files_full": listed_f, "stage_table_source": source, "total_ps": r6(T_eq),
           "n_frames_solute": NS, "n_frames_full": NFe, "files_written": written}
    if NS:
        out["imaging_solute"] = {k: acc[k] for k in ("frames_with_box", "frames_without_box",
                                                     "frames_zinc_moved", "frames_ligand_moved")}
    out["state_of_start"] = start_state
    out["n2_sign_of_start"] = start_sign
    out["stages"] = stage_list
    if NS:
        stt = list(cols["state"])
        sg = [int(x) for x in cols["n2_sign"]]
        out.update({"count_C": stt.count("C"), "count_R": stt.count("R"), "count_I": stt.count("I"),
                    "count_n2_positive": sg.count(1), "count_n2_negative": sg.count(-1)})
        i1 = next((i for i in range(NS) if stt[i] != start_state), None)
        i2 = next((i for i in range(NS) if stt[i] == opposite), None)
        i3 = next((i for i in range(NS) if sg[i] != start_sign), None)
        out["first_frame_not_in_state_of_start"] = frame_entry(i1) if i1 is not None else None
        out["first_frame_in_opposite_state"] = frame_entry(i2) if i2 is not None else None
        out["first_frame_n2_sign_not_of_start"] = frame_entry(i3) if i3 is not None else None
        out["state_last_frame"] = stt[-1]
        out["n2_sign_last_frame"] = sg[-1]
        out["time_last_solute_frame_ns"] = r6(lead[-1][4])
    else:
        out.update({"count_C": 0, "count_R": 0, "count_I": 0, "count_n2_positive": 0, "count_n2_negative": 0,
                    "first_frame_not_in_state_of_start": None, "first_frame_in_opposite_state": None,
                    "first_frame_n2_sign_not_of_start": None, "state_last_frame": None,
                    "n2_sign_last_frame": None, "time_last_solute_frame_ns": None})
    out["_aux"] = eaux
    return out


# ----------------------------------------------------------------------------------------------
# campaign (section 12)
# ----------------------------------------------------------------------------------------------
def find_run_folders(top):
    found = []
    for root, dirs, files in os.walk(top):
        dirs.sort()
        dirs[:] = [d for d in dirs if not d.startswith("analysis")]
        if any(SEG_RE.match(f) and f.startswith("solute_") for f in files):
            found.append(root)
    return sorted(found)


def campaign_rules(rows):
    """rows: table rows of the production runs of one zinc model."""
    by = {s: [r for r in rows if r["start"] == s] for s in "XNR"}
    n_left = sum(1 for r in by["X"] if r["rule1_left_crystal"])
    rule1_failed = n_left >= 2
    conv = {s: sum(1 for r in by[s] if r["verdict"] == "CONVERTED") for s in "XNR"}
    stay = {s: sum(1 for r in by[s] if r["verdict"] == "STAYED") for s in "XNR"}
    mixed = {s: sum(1 for r in by[s] if r["verdict"] == "MIXED") for s in "XNR"}
    if rule1_failed:
        outcome = "validation_failed"
    elif conv["R"] >= 2 and all(r["verdict"] == "STAYED" for r in by["X"] + by["N"]):
        outcome = "rotated_converts"
    elif all(r["verdict"] == "STAYED" for r in rows):
        outcome = "both_persist"
    else:
        outcome = "mixed"
    complete = {s: sum(1 for r in by[s] if r["complete"]) for s in "XNR"}
    nine = all(len(by[s]) == 3 and complete[s] == 3 for s in "XNR") and len(rows) == 9
    return {"n_runs": len(rows), "n_complete": sum(1 for r in rows if r["complete"]),
            "runs_per_start": {s: len(by[s]) for s in "XNR"}, "all_nine_complete": nine,
            "final": bool(all(complete[s] >= 3 for s in "XNR")),
            "rule1": {"n_X_runs": len(by["X"]), "n_X_left_crystal": n_left, "failed": rule1_failed},
            "n_converted": conv, "n_stayed": stay, "n_mixed": mixed, "outcome": outcome,
            "rule4_sensitivity_run_required": outcome == "rotated_converts"}


OUTCOME_TEXT = {
    "validation_failed": "rule 1 fails: the crystal-start runs leave the crystal orientation; the simulation "
                         "supports nothing and is reported as a failed validation",
    "rotated_converts": "the rotated starts convert, the crystal and near-native starts stay: the rotated "
                        "orientation is not stable outside the scoring function",
    "both_persist": "no run converts: both orientations persisted; the simulation does not decide",
    "mixed": "mixed: the runs are reported one by one; no average over runs that disagree"}


def run_campaign(top, args, summaries=None, out_dir=None):
    """Analyse the production runs below `top` (or take summaries already made: {run_dir: summary})
    and evaluate the rules.  Returns the campaign record."""
    top = os.path.abspath(top)
    out_dir = out_dir or (os.path.abspath(args.out) if getattr(args, "out", None)
                          else os.path.join(top, "analysis_campaign_recheck"))
    rows, skipped, failed, runs = [], [], [], []
    for d in find_run_folders(top):
        rel = os.path.relpath(d, top)
        name = os.path.basename(d)
        if not CAMPAIGN_NAME_RE.match(name):
            skipped.append({"folder": rel, "reason": "the folder name is not a production run (<start><replica>, "
                                                     "replica 1, 2, ...)"})
            continue
        try:
            if summaries is not None and d in summaries:
                summ, odir = summaries[d]
            else:
                a = argparse.Namespace(**vars(args))
                a.out = None
                a.start = a.model = a.planned_ns = a.system = None
                summ, odir = analyse_run(d, a, quiet=True)
        except Stop as err:
            failed.append({"folder": rel, "error": str(err)})
            continue
        if summ["run"]["role"] == "validation":
            skipped.append({"folder": rel, "reason": "the run has the role validation by its manifest"})
            continue
        spath = os.path.join(odir, "summary.json")
        runs.append({"table_row": summ["table_row"], "summary_file": spath, "summary_sha256": sha256_file(spath),
                     "folder": rel})
    order = {"X": 0, "N": 1, "R": 2}

    def key(r):
        t = r["table_row"]
        rep = t["replica"] if isinstance(t["replica"], int) else 10 ** 9
        return (str(t["model"]), order[t["start"]], rep, r["folder"])
    runs.sort(key=key)
    rows = [r["table_row"] for r in runs]
    models = []
    for r in rows:
        if r["model"] not in models:
            models.append(r["model"])
    rules = {m: campaign_rules([r for r in rows if r["model"] == m]) for m in models}
    rec = {"spec": "md/analysis/SPEC.md version " + SPEC_VERSION,
           "script": {"name": PROGRAM, "sha256": sha256_file(os.path.abspath(__file__))},
           "software": {"python": platform.python_version(), "numpy": np.__version__},
           "campaign_dir": top, "rules_per_model": rules,
           "outcome_text": {m: OUTCOME_TEXT[rules[m]["outcome"]] + ("" if rules[m]["final"] else " (not final)")
                            for m in models},
           "runs": runs, "failed": failed, "skipped": skipped}
    os.makedirs(out_dir, exist_ok=True)
    write_json(os.path.join(out_dir, "campaign.json"), rec)
    with open(os.path.join(out_dir, "campaign_table.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(TABLE_ROW_KEYS)
        for r in rows:
            w.writerow(["" if r[k] is None else r[k] for k in TABLE_ROW_KEYS])
    return rec, out_dir


# ----------------------------------------------------------------------------------------------
# comparison of two implementations (section 14)
# ----------------------------------------------------------------------------------------------
TOL = {"dist": 0.001, "ang": 0.01, "tors": 0.01, "vol": 0.001, "time_ns": 0.000001, "time_ps": 0.001,
       "number": 0.000001}
EPS = 1e-9
SKIP_KEYS = {"note", "file", "run_dir", "system_dir", "summary_file", "summary_sha256", "files_written",
             "files_solute", "files_full"}
# texts whose wording the specification does not fix: compared, listed apart, never a disagreement
LABEL_KEYS = {"order", "kind", "trajectory", "label", "stage_table_source"}
LABEL_PATHS = {"table_row/run"}
# the number of warnings is a key of table_row and is compared (identical); the wording is not (section 14)
WARNING_PATHS = {"table_row/n_warnings"}
THRESHOLDS = {
    "hb_N2H_T199OG1": [("d_N2_T199OG1_A", HB_DIST_A, 0.001), ("ang_N2_HN2_T199OG1_deg", HB_ANGLE_DEG, 0.01)],
    "hb_T199NH_Osulfonyl": [("d_T199N_O1_A", HB_DIST_A, 0.001), ("d_T199N_O2_A", HB_DIST_A, 0.001),
                            ("ang_T199N_H_O1_deg", HB_ANGLE_DEG, 0.01), ("ang_T199N_H_O2_deg", HB_ANGLE_DEG, 0.01)],
    "hb_T199OG1H_E106OE": [("d_T199OG1_E106OE1_A", HB_DIST_A, 0.001), ("d_T199OG1_E106OE2_A", HB_DIST_A, 0.001),
                           ("ang_T199OG1_HG1_E106OE1_deg", HB_ANGLE_DEG, 0.01),
                           ("ang_T199OG1_HG1_E106OE2_deg", HB_ANGLE_DEG, 0.01)],
    "fifth_sulfonylO": [("d_Zn_Osulfonyl_min_A", FIFTH_A, 0.001)],
    "fifth_E106OE1": [("d_Zn_E106OE1_A", FIFTH_A, 0.001)],
    "fifth_E106OE2": [("d_Zn_E106OE2_A", FIFTH_A, 0.001)],
    "fifth_E106OE": [("d_Zn_E106OE_min_A", FIFTH_A, 0.001)],
    "state": [("d1_N1_F131CZ_A", MID_D1_A, 0.001), ("d2_N1_N67OD1_A", MID_D2_A, 0.001)],
    "n2_sign": [("n2_volume", 0.0, 0.001)],
    "fifth_water": [("d_Zn_Owater_min_A", FIFTH_A, 0.001)],
    "fifth_proteinO": [("d_Zn_Oprotein_min_A", FIFTH_A, 0.001)],
    "fifth_any": [("d_Zn_O_any_min_A", FIFTH_A, 0.001)],
}
FOLLOWS = {"state_assigned": "state", "n2_sign_assigned": "n2_sign"}
# Which threshold frames can explain a difference under a key of summary.json (section 14: "counts may
# differ by at most their number, and a verdict that changes with them is reported as it is with both
# readings").  A difference in a block is listed as following from threshold frames only if a column
# named here has such frames in this run; None = any column (the table row gathers all blocks).
GATE_FLAGS = {"hb_N2H_T199OG1", "fifth_sulfonylO", "fifth_any", "fifth_water", "fifth_proteinO"}
THRESHOLD_SOURCES = {
    "lig_rmsd_labelling_counts": {"lig_rmsd_labelling"},
    "orientation": {"state", "state_assigned"},
    "asn67_ND2": {"asn67_ND2_nearer_N1", "state"},
    "n2_configuration": {"n2_sign", "n2_sign_assigned", "state"},
    "hbond_first_passage": {"state", "state_assigned", "hb_N2H_T199OG1"},
    "start_of_production": {"state", "n2_sign"},
    "equilibration": {"state", "n2_sign", "n_water_ring_place", "n_water_ring_wide_3p5A"},
    "validation": {"hb_N2H_T199OG1"},
    "gate_GV": GATE_FLAGS | {"fifth_E106OE1", "fifth_E106OE2"},
    "table_row": None,
    "stats": set(), "run": set(), "water_ring_place": set(), "rmsf": set(),
}
# keys under which a difference, even one that section 14 allows, must be read by a person: the
# comparison then ends with exit status 3
VERDICT_KEYS = {"verdict", "outcome", "passed", "rule1_left_crystal", "gate_GV_passed", "complete", "role",
                "provisional", "failed", "final", "rule4_sensitivity_run_required"}


def csv_kind(col):
    if col == "time_ns":
        return "time_ns"
    if col.startswith("tors_"):
        return "tors"
    if col.endswith("_deg"):
        return "ang"
    if col.endswith("_A"):
        return "dist"
    if col == "n2_volume":
        return "vol"
    if col in ("state", "state_assigned", "residue_name"):
        return "text"
    return "int"


def read_csv(path):
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh))
    if not rows:
        return [], []
    return rows[0], rows[1:]


def compare_csv(name, p_path, r_path, aux):
    """Line by line.  p = primary (analyse_md.py), r = this program."""
    res = {"file": name, "primary": p_path, "recheck": r_path, "present_primary": os.path.isfile(p_path),
           "present_recheck": os.path.isfile(r_path), "columns": {}, "disagreements": [], "threshold_frames": [],
           "n_values": 0, "n_disagree": 0, "n_threshold": 0, "notes": []}
    if not res["present_primary"] and not res["present_recheck"]:
        res["agree"] = True
        res["notes"].append("written by neither program")
        return res
    if res["present_primary"] != res["present_recheck"]:
        res["agree"] = False
        res["n_disagree"] = 1
        res["disagreements"].append({"what": "file written by one program only"})
        return res
    hp, rp = read_csv(p_path)
    hr, rr = read_csv(r_path)
    res["lines_primary"], res["lines_recheck"] = len(rp), len(rr)
    res["columns_only_primary"] = [c for c in hp if c not in hr]
    res["columns_only_recheck"] = [c for c in hr if c not in hp]
    res["column_order_equal"] = hp == hr
    common = [c for c in hr if c in hp]
    if len(rp) != len(rr):
        res["n_disagree"] += 1
        res["disagreements"].append({"what": "number of lines", "primary": len(rp), "recheck": len(rr)})
    if res["columns_only_primary"] or res["columns_only_recheck"] or not res["column_order_equal"]:
        res["n_disagree"] += 1
        res["disagreements"].append({"what": "columns", "only_primary": res["columns_only_primary"],
                                     "only_recheck": res["columns_only_recheck"],
                                     "same_order": res["column_order_equal"]})
    ip = {c: hp.index(c) for c in common}
    ir = {c: hr.index(c) for c in common}
    n = min(len(rp), len(rr))
    excused_cols = {}
    tie_lab = set(aux.get("labelling_tie_frames", [])) if name == "timeseries_solute.csv" else set()
    pre = {"timeseries_full.csv": "full_frames_", "timeseries_equil_full.csv": "equil_full_frames_"}.get(name)
    if pre:
        edge_water = {"n_water_N1_3p5A": set(aux.get(pre + "water_N1_at_3p5A", [])),
                      "n_water_ring_wide_3p5A": set(aux.get(pre + "water_ring_wide_at_3p5A", [])),
                      "n_water_ring_place": set(aux.get(pre + "water_ring_place_at_cutoff", []))}
    else:
        edge_water = {}

    def val(row, idx, c):
        return row[idx[c]] if idx[c] < len(row) else ""

    def near(i, quantity, thr, tol):
        if quantity not in ip:
            return False
        for row, idx in ((rp[i], ip), (rr[i], ir)):
            try:
                if abs(float(val(row, idx, quantity)) - thr) <= tol + 1e-4:
                    return True
            except ValueError:
                pass
        return False
    for c in common:
        kind = csv_kind(c)
        st = {"kind": kind, "n": n, "n_differ": 0, "n_threshold": 0}
        if kind in TOL:
            st["tolerance"] = TOL[kind]
            st["max_abs_diff"] = 0.0
        for i in range(n):
            a, b = val(rp[i], ip, c), val(rr[i], ir, c)
            res["n_values"] += 1
            if kind in TOL:
                try:
                    fa, fb = float(a), float(b)
                except ValueError:
                    same = a == b
                    diff = None
                else:
                    diff = abs(fa - fb)
                    if kind == "tors":
                        diff = abs((fa - fb + 180.0) % 360.0 - 180.0)
                    if math.isnan(fa) and math.isnan(fb):
                        diff = 0.0
                    same = diff <= TOL[kind] + EPS
                    if diff == diff:
                        st["max_abs_diff"] = max(st["max_abs_diff"], diff)
                if not same:
                    st["n_differ"] += 1
                    if len(res["disagreements"]) < 40:
                        res["disagreements"].append({"line": i, "column": c, "primary": a, "recheck": b, "diff": diff,
                                                     "tolerance": TOL[kind]})
                continue
            if kind == "int":
                try:
                    same = int(float(a)) == int(float(b))
                except ValueError:
                    same = a == b
            else:
                same = a == b
            if same:
                continue
            excuse = None
            if c in THRESHOLDS and any(near(i, q, t, tol) for q, t, tol in THRESHOLDS[c]):
                excuse = "value within the tolerance of its threshold"
            elif c in FOLLOWS and any(x[0] <= i for x in excused_cols.get(FOLLOWS[c], [])):
                excuse = "follows from a threshold frame of " + FOLLOWS[c]
            elif c == "asn67_ND2_nearer_N1" and "d_N1_N67ND2_A" in ip and "d2_N1_N67OD1_A" in ip:
                try:
                    if abs(float(val(rr[i], ir, "d_N1_N67ND2_A")) - float(val(rr[i], ir, "d2_N1_N67OD1_A"))) < 0.001 + 1e-4:
                        excuse = "the two distances differ by less than 0.001 A"
                except ValueError:
                    pass
            elif c == "lig_rmsd_labelling" and i in tie_lab:
                excuse = "two labellings give the same RMSD within 0.001 A"
            elif c in edge_water and i in edge_water[c]:
                excuse = "a water oxygen within 0.001 A of the cut-off of the count"
            elif c in ("water_O_index_nearest_Zn", "protein_O_index_nearest_Zn"):
                q = "d_Zn_Owater_min_A" if c.startswith("water") else "d_Zn_Oprotein_min_A"
                try:
                    if abs(float(val(rp[i], ip, q)) - float(val(rr[i], ir, q))) <= 0.001:
                        excuse = "two atoms equally near within 0.001 A (not verified atom by atom)"
                except (ValueError, KeyError):
                    pass
            if excuse:
                st["n_threshold"] += 1
                try:
                    dv = abs(int(float(a)) - int(float(b)))
                except ValueError:
                    dv = None
                if dv is not None:
                    st["threshold_abs_sum"] = st.get("threshold_abs_sum", 0) + dv
                    st["threshold_sq_sum"] = st.get("threshold_sq_sum", 0) + dv * dv
                    st["threshold_abs_max"] = max(st.get("threshold_abs_max", 0), dv)
                excused_cols.setdefault(c, []).append((i, excuse))
                if len(res["threshold_frames"]) < 200:
                    res["threshold_frames"].append({"line": i, "column": c, "primary": a, "recheck": b, "why": excuse})
            else:
                st["n_differ"] += 1
                if len(res["disagreements"]) < 40:
                    res["disagreements"].append({"line": i, "column": c, "primary": a, "recheck": b})
        res["columns"][c] = st
        res["n_disagree"] += st["n_differ"]
        res["n_threshold"] += st["n_threshold"]
    res["agree"] = res["n_disagree"] == 0
    return res


def key_kind(path):
    """Tolerance class of a number in summary.json from its place."""
    last = next((p for p in reversed(path) if isinstance(p, str)), "")
    if len(path) >= 2 and path[0] == "stats":
        col = path[1]
        if last in ("n", "n_last_half", "block_frames"):
            return "exact"
        k = column_kind(col)
        return {"count": "number"}.get(k, k)
    if last.endswith("_ns"):
        return "time_ns"
    if last.endswith("_ps"):
        return "time_ps"
    if last.endswith("_A") or last in ("window_A",):
        return "dist"
    if last.endswith("_deg"):
        return "ang"
    if last.startswith("n2_volume"):
        return "vol"
    if "fraction" in last or last.startswith("frac_"):
        return "fraction"
    return "number"


def compare_tree(p, r, path, out, ctx):
    """Recursive comparison of two JSON values; p = primary, r = recheck."""
    last = path[-1] if path else ""
    if isinstance(last, str) and last in SKIP_KEYS:
        out["n_skipped"] += 1
        return
    spath = "/".join(map(str, path))
    if (isinstance(last, str) and last in LABEL_KEYS) or spath in LABEL_PATHS:
        out["n_labels"] = out.get("n_labels", 0) + 1
        if p != r:
            out.setdefault("labels_differ", []).append({"key": spath, "primary": p, "recheck": r})
        return
    if spath in WARNING_PATHS:
        out["n_compared"] += 1
        if p != r:
            out["disagreements"].append({"key": spath, "primary": p, "recheck": r, "diff": None, "tolerance": 0,
                                         "kind": "number of warnings (the texts are listed below, their wording is not compared)"})
        return
    if isinstance(p, dict) and isinstance(r, dict):
        for k in p:
            if k not in r:
                if k in SKIP_KEYS:
                    continue
                out["only_primary"].append("/".join(map(str, path + [k])))
        for k in r:
            if k not in p:
                if k in SKIP_KEYS or k.startswith("_"):
                    continue
                out["only_recheck"].append("/".join(map(str, path + [k])))
        for k in r:
            if k in p:
                compare_tree(p[k], r[k], path + [k], out, ctx)
        return
    if isinstance(p, list) and isinstance(r, list):
        if len(p) != len(r):
            record(out, path, p if len(p) < 8 else "list of %d" % len(p), r if len(r) < 8 else "list of %d" % len(r),
                   None, None, "length of a list", ctx)
            return
        for i, (a, b) in enumerate(zip(p, r)):
            compare_tree(a, b, path + [i], out, ctx)
        return
    out["n_compared"] += 1
    if p is None or r is None or isinstance(p, (bool, str)) or isinstance(r, (bool, str)):
        if p != r or type(p) is not type(r):
            if isinstance(p, bool) != isinstance(r, bool) and p == r:
                return
            record(out, path, p, r, None, None, "identical", ctx)
        return
    if isinstance(p, (int, float)) and isinstance(r, (int, float)):
        kind = key_kind(path)
        if spath in EXACT_PATHS:
            if abs(p - r) > 1e-12:
                record(out, path, p, r, abs(p - r), 0, "identical", ctx)
            return
        if kind == "exact" or (isinstance(p, int) and isinstance(r, int) and kind in ("number", "exact")):
            if p != r:
                record(out, path, p, r, abs(p - r), 0, "identical", ctx)
            return
        tol = {"fraction": 0.000001}.get(kind, TOL.get(kind, 0.000001))
        diff = abs(p - r)
        if kind == "tors" and last in ("mean", "mean_last_half") or (kind == "tors" and isinstance(last, int)):
            diff = abs((p - r + 180.0) % 360.0 - 180.0)
        out["max_diff"][kind] = max(out["max_diff"].get(kind, 0.0), diff)
        if diff > tol + EPS:
            record(out, path, p, r, diff, tol, kind, ctx)
        return
    record(out, path, p, r, None, None, "type", ctx)


def record(out, path, p, r, diff, tol, kind, ctx):
    spath = "/".join(map(str, path))
    item = {"key": spath, "primary": p, "recheck": r, "diff": diff, "tolerance": tol, "kind": kind}
    # a part of gate G-V or rule 1 whose mean lies within the tolerance of its limit
    why = None
    par = ctx.get("parent_of", {}).get(spath)
    if isinstance(par, dict) and path and path[-1] == "passed":
        for mkey, wkey, t in (("mean_A", "window_A", 0.001), ("mean_deg", "window_deg", 0.01)):
            if is_number(par.get(mkey)) and isinstance(par.get(wkey), list):
                if any(abs(par[mkey] - b) <= t for b in par[wkey]):
                    why = "mean within the tolerance of a limit of gate G-V"
        if is_number(par.get("mean_A")) and is_number(par.get("limit_max_A")) and abs(par["mean_A"] - par["limit_max_A"]) <= 0.001:
            why = "mean within the tolerance of a limit of gate G-V"
    if isinstance(par, dict) and path and path[-1] == "rule1_left_crystal":
        m_ = par.get("lig_rmsd_mean_last_half_A")       # the same key in validation and in table_row
        if is_number(m_) and abs(m_ - RULE1_A) <= 0.001:
            why = "mean of the last half within the tolerance of the limit of rule 1 (2.0 A)"
    if why is None and path and path[0] == "gate_GV" and path[-1] in ("passed", "verdict") and ctx.get("gate_part_at_limit"):
        why = "follows from a part of gate G-V whose mean lies within the tolerance of its limit"
    if why is None and path and path[0] == "table_row" and path[-1] == "gate_GV_passed" and ctx.get("gate_part_at_limit"):
        why = "follows from a part of gate G-V whose mean lies within the tolerance of its limit"
    nt = ctx.get("n_threshold", 0)
    # the statistics of a water count of the full trajectory (and its mean in the table row) when the two
    # programs class a water oxygen at the cut-off differently in some frames (section 14: such frames are
    # listed and counts may differ by at most their number): the difference must stay within what these
    # frames can make, else it is a disagreement
    ccol = None
    if len(path) >= 3 and path[0] == "stats" and path[1] in FULL_COUNTS:
        ccol, ckey = path[1], path[2]
    elif len(path) == 2 and path[0] == "table_row" and path[1].endswith("_mean") and path[1][:-5] in FULL_COUNTS:
        ccol, ckey = path[1][:-5], "mean"
    th = (ctx.get("count_threshold") or {}).get(ccol) if ccol else None
    if th and is_number(p) and is_number(r) and why is None:
        stt = ((ctx.get("recheck_tree") or {}).get("stats") or {}).get(ccol) or {}
        n_all, n_half, blk = stt.get("n") or 0, stt.get("n_last_half") or 0, stt.get("block_frames") or 0
        bound = None
        if ckey == "mean" and n_all:
            bound = th["abs_sum"] / n_all
        elif ckey == "mean_last_half" and n_half:
            bound = th["abs_sum"] / n_half
        elif ckey == "sd" and n_all:
            bound = math.sqrt(th["sq_sum"] / n_all)
        elif ckey == "sd_last_half" and n_half:
            bound = math.sqrt(th["sq_sum"] / n_half)
        elif ckey in ("min", "max"):
            bound = th["abs_max"]
        elif ckey in ("block_means", "block_sem") and blk:
            bound = th["abs_sum"] / blk
        if bound is not None and abs(p - r) <= bound + 1e-6:
            item["why"] = ("follows from the %d threshold frame(s) of %s (a water oxygen within 0.001 A of the "
                           "cut-off): at most %.6g; both readings are given" % (th["n"], ccol, bound))
            out["threshold"].append(item)
            return
    if len(path) == 5 and path[0] == "equilibration" and path[1] == "stages" and path[3] == "driver_hydration_check":
        # section 14: count_agrees may differ only through a water oxygen within 0.001 A of the
        # cut-off in the frame compared, nearest_agrees only if the difference from the driver's
        # distance lies within 0.0001 A of 0.002 A; every other key of the block must be identical
        tag = None
        try:
            tag = ctx["recheck_tree"]["equilibration"]["stages"][path[2]]["tag"]
        except (KeyError, IndexError, TypeError):
            pass
        hint = (ctx.get("driver_check") or {}).get(tag) or {}
        if path[4] == "count_agrees" and hint.get("water_oxygen_at_cutoff"):
            item["why"] = "a water oxygen within 0.001 A of the cut-off in the frame compared with the driver"
        elif path[4] == "nearest_agrees" and is_number(hint.get("abs_difference_from_driver_A")) and \
                abs(hint["abs_difference_from_driver_A"] - DRIVER_NEAREST_TOL_A) <= 0.0001:
            item["why"] = "the difference from the driver's distance lies within 0.0001 A of 0.002 A"
        if "why" in item:
            out["threshold"].append(item)
        else:
            out["disagreements"].append(item)
        return
    thr_cols = ctx.get("threshold_columns") or {}
    block = path[0] if path else None
    if block == "fractions" and len(path) >= 2:
        srcs = {path[1]}
    else:
        srcs = THRESHOLD_SOURCES.get(block, set())
    nt_rel = nt if srcs is None else sum(n_ for c_, n_ in thr_cols.items() if c_ in srcs)
    if why is None and nt_rel > 0 and kind in ("identical", "fraction", "length of a list"):
        if is_number(p) and is_number(r) and kind == "identical" and abs(p - r) > nt_rel:
            why = None
        else:
            why = "may follow from the %d threshold frame(s) of this run in %s; both readings are given" % (
                nt_rel, "any column" if srcs is None else ", ".join(sorted(c_ for c_ in srcs if thr_cols.get(c_))))
    if why:
        item["why"] = why
        out["threshold"].append(item)
    else:
        out["disagreements"].append(item)


def parents(tree, path=(), acc=None):
    acc = {} if acc is None else acc
    if isinstance(tree, dict):
        for k, v in tree.items():
            acc["/".join(map(str, path + (k,)))] = tree
            parents(v, path + (k,), acc)
    elif isinstance(tree, list):
        for i, v in enumerate(tree):
            parents(v, path + (i,), acc)
    return acc


COMPARED_BLOCKS = ["stats", "fractions", "lig_rmsd_labelling_counts", "orientation", "asn67_ND2",
                   "n2_configuration", "hbond_first_passage", "water_ring_place", "start_of_production",
                   "equilibration", "validation", "gate_GV", "rmsf", "table_row"]
# numbers that section 14 wants identical although their key ends like a distance
EXACT_PATHS = {"water_ring_place/water_cutoff_A"}
RUN_KEYS = ["n_frames_solute", "n_frames_full", "seen_ns", "planned_ns", "complete", "role"]
CSV_FILES = ["timeseries_solute.csv", "timeseries_full.csv", "rmsf_ca.csv", "timeseries_equil_solute.csv",
             "timeseries_equil_full.csv"]


def strip_equilibration(e):
    """Section 14: the equilibration block without the file and header entries of its stages."""
    if not isinstance(e, dict):
        return e
    e = dict(e)
    if isinstance(e.get("stages"), list):
        e["stages"] = [{k: v for k, v in st.items() if k not in ("solute", "full")} if isinstance(st, dict) else st
                       for st in e["stages"]]
    return e


def compare_run(primary_dir, recheck_dir):
    """Compare the output of analyse_md.py (primary_dir) with that of this program (recheck_dir)."""
    res = {"primary_dir": primary_dir, "recheck_dir": recheck_dir, "files": {}, "summary": None}
    ps, rs = os.path.join(primary_dir, "summary.json"), os.path.join(recheck_dir, "summary.json")
    if not os.path.isfile(ps):
        res["agree"] = None
        res["error"] = "no summary.json of the primary implementation in %s" % primary_dir
        return res
    with open(ps) as fh:
        P = json.load(fh)
    with open(rs) as fh:
        Rr = json.load(fh)
    if spec_version_of(P) != SPEC_VERSION:
        res["agree"] = None
        res["error"] = ("the output of the first program in %s was written after %s; this program follows version %s "
                        "(run analyse_md.py again, or give --run-primary)" % (primary_dir, P.get("spec"), SPEC_VERSION))
        return res
    aux = {}
    ap = os.path.join(recheck_dir, "recheck_aux.json")
    if os.path.isfile(ap):
        with open(ap) as fh:
            aux = json.load(fh)
    n_thr = 0
    for name in CSV_FILES:
        c = compare_csv(name, os.path.join(primary_dir, name), os.path.join(recheck_dir, name), aux)
        res["files"][name] = c
        n_thr += c["n_threshold"]
    out = {"n_compared": 0, "n_skipped": 0, "disagreements": [], "threshold": [], "only_primary": [],
           "only_recheck": [], "max_diff": {}}
    pt = {b: P.get(b) for b in COMPARED_BLOCKS}
    rt = {b: Rr.get(b) for b in COMPARED_BLOCKS}
    pt["equilibration"], rt["equilibration"] = strip_equilibration(pt["equilibration"]), strip_equilibration(rt["equilibration"])
    pt["run"] = {k: mget(P, "run", k) for k in RUN_KEYS}
    rt["run"] = {k: mget(Rr, "run", k) for k in RUN_KEYS}
    count_thr = {}
    for col, st in (res["files"].get("timeseries_full.csv") or {}).get("columns", {}).items():
        if col in FULL_COUNTS and st.get("n_threshold"):
            count_thr[col] = {"n": st["n_threshold"], "abs_sum": st.get("threshold_abs_sum", 0),
                              "sq_sum": st.get("threshold_sq_sum", 0), "abs_max": st.get("threshold_abs_max", 0)}
    thr_cols = {}
    for f_ in res["files"].values():
        for col, st in f_.get("columns", {}).items():
            if st.get("n_threshold"):
                thr_cols[col] = thr_cols.get(col, 0) + st["n_threshold"]
    # a part of gate G-V whose mean lies within the tolerance of its limit (either program's mean)
    at_limit = []
    for tree in (pt.get("gate_GV"), rt.get("gate_GV")):
        for key, par in parents(tree or {}).items():
            if not key.endswith("/passed") or not isinstance(par, dict):
                continue
            for mkey, wkey, t in (("mean_A", "window_A", 0.001), ("mean_deg", "window_deg", 0.01)):
                if is_number(par.get(mkey)) and isinstance(par.get(wkey), list) and \
                        any(is_number(b) and abs(par[mkey] - b) <= t for b in par[wkey]):
                    at_limit.append(key)
            if is_number(par.get("mean_A")) and is_number(par.get("limit_max_A")) and \
                    abs(par["mean_A"] - par["limit_max_A"]) <= 0.001:
                at_limit.append(key)
    ctx = {"n_threshold": n_thr, "parent_of": parents(rt), "recheck_tree": rt,
           "driver_check": aux.get("equil_driver_check", {}), "count_threshold": count_thr,
           "threshold_columns": thr_cols, "gate_part_at_limit": sorted(set(at_limit))}
    compare_tree(pt, rt, [], out, ctx)
    out["spec_primary"] = P.get("spec")
    out["script_primary"] = P.get("script")
    here = primary_script()
    out["script_primary_is_the_file_beside_this_program"] = (
        sha256_file(here) == mget(P, "script", "sha256") if os.path.isfile(here) else None)
    out["warnings_primary"] = P.get("warnings")
    out["warnings_recheck"] = Rr.get("warnings")
    out["n_warnings_primary"] = len(P.get("warnings") or [])
    out["n_warnings_recheck"] = len(Rr.get("warnings") or [])
    out["labels"] = {"primary": mget(P, "run", "label"), "recheck": mget(Rr, "run", "label")}
    res["summary"] = out
    res["n_disagree"] = len(out["disagreements"]) + sum(c["n_disagree"] for c in res["files"].values())
    res["n_threshold"] = len(out["threshold"]) + n_thr
    res["n_keys_one_side"] = len(out["only_primary"]) + len(out["only_recheck"])
    res["n_values"] = out["n_compared"] + sum(c["n_values"] for c in res["files"].values())
    res["warning_count_equal"] = out["n_warnings_primary"] == out["n_warnings_recheck"]
    res["threshold_columns"] = thr_cols
    res["gate_parts_at_a_limit"] = ctx["gate_part_at_limit"]
    res["verdict_depends_on_threshold_frames"] = any(
        d["key"].split("/")[-1] in VERDICT_KEYS for d in out["threshold"])
    res["keys_that_depend_on_threshold_frames"] = [d["key"] for d in out["threshold"]
                                                   if d["key"].split("/")[-1] in VERDICT_KEYS]
    res["agree"] = res["n_disagree"] == 0 and res["n_keys_one_side"] == 0
    return res


def skipped_names(lst, top=None):
    """The folders a campaign left out, as a set of paths relative to the campaign folder `top`,
    whatever the form of the entries (a path, or an object that holds one).  A path that is absolute
    is taken relative to `top`; if it does not lie below `top`, its last two parts are taken."""
    out = set()
    top = os.path.realpath(top) if isinstance(top, str) and top else None
    for e in lst or []:
        if isinstance(e, str):
            vals = [e]
        elif isinstance(e, dict):
            vals = [e[k] for k in ("folder", "run_dir", "dir", "path", "run", "label", "name") if isinstance(e.get(k), str)]
            vals += [v for v in e.values() if isinstance(v, str) and v not in vals]
        else:
            vals = []
        for v in vals:
            parts = [x for x in re.split(r"[\\/]+", v.strip()) if x]
            if not parts or " " in v:
                continue
            if os.path.isabs(v):
                rel = os.path.relpath(os.path.realpath(v), top) if top else None
                out.add(rel.replace(os.sep, "/") if rel and not rel.startswith("..") else "/".join(parts[-2:]))
            else:
                out.add("/".join(parts))
            break
    return out


def compare_campaign(primary_json, rec):
    res = {"primary": primary_json, "present_primary": os.path.isfile(primary_json)}
    if not res["present_primary"]:
        res["agree"] = None
        res["error"] = "no campaign.json of the primary implementation: %s" % primary_json
        return res
    with open(primary_json) as fh:
        P = json.load(fh)
    out = {"n_compared": 0, "n_skipped": 0, "disagreements": [], "threshold": [], "only_primary": [],
           "only_recheck": [], "max_diff": {}}
    compare_tree(P.get("rules_per_model"), rec["rules_per_model"], ["rules_per_model"], out, {})
    sp = skipped_names(P.get("skipped"), P.get("campaign_dir") or rec.get("campaign_dir"))
    sr = skipped_names(rec["skipped"], rec.get("campaign_dir"))
    out["skipped_primary"], out["skipped_recheck"] = sorted(sp), sorted(sr)
    out["n_compared"] += 1
    if sp != sr:
        out["disagreements"].append({"key": "skipped", "primary": sorted(sp), "recheck": sorted(sr),
                                     "kind": "identical"})
    pr = [mget(x, "table_row", "run") for x in P.get("runs") or []]
    rr = [mget(x, "table_row", "run") for x in rec["runs"]]
    out["runs_primary"], out["runs_recheck"] = pr, rr
    out["n_compared"] += 1
    if [str(x).split("/")[-1] for x in pr] != [str(x).split("/")[-1] for x in rr]:
        out["disagreements"].append({"key": "runs (order)", "primary": pr, "recheck": rr, "kind": "identical"})
    res.update(out)
    res["n_disagree"] = len(out["disagreements"])
    res["agree"] = res["n_disagree"] == 0 and not out["only_primary"] and not out["only_recheck"]
    return res


def report_markdown(rep):
    L = []
    L.append("# Second analysis against the first: %s" % rep["what"])
    L.append("")
    L.append("Written by `md/analysis/recheck_md.py` (sha256 `%s`); first program `analyse_md.py` (sha256 `%s`). "
             "Specification: %s. Tolerances: section 14 (distances 0.001 Å, angles 0.01°, signed volume 0.001, "
             "times 0.000001 ns, counts, flags, states and verdicts identical)."
             % (rep["script_sha256"][:16] + "…", (rep["primary_script_sha256"] or "not found")[:16] + "…", rep["spec"]))
    L.append("")
    L.append("**Result: %s.**" % rep["result_text"])
    L.append("")
    for v in rep.get("verdicts_that_differ_through_threshold_frames") or []:
        L.append("- **To be read**, run %s: %s differ between the two programs through frames at a threshold "
                 "or a mean at a limit (both readings below)." % (v["run"], ", ".join("`%s`" % k for k in v["keys"])))
    if rep.get("verdicts_that_differ_through_threshold_frames"):
        L.append("")
    L.append("| Run | Solute / full frames | Values compared | Disagreements | Threshold frames and what follows | "
             "Keys in one program only | Warnings (first / second) | Agree |")
    L.append("|---|---|---|---|---|---|---|---|")
    for r in rep["runs"]:
        if r.get("error"):
            L.append("| %s | | | | | | | **not compared**: %s |" % (r["label"], r["error"]))
            continue
        c = r["comparison"]
        L.append("| %s | %s / %s | %d | %d | %d | %d | %d / %d | %s |"
                 % (r["label"], r["n_frames_solute"], r["n_frames_full"], c["n_values"], c["n_disagree"],
                    c["n_threshold"], c["n_keys_one_side"], c["summary"]["n_warnings_primary"],
                    c["summary"]["n_warnings_recheck"], "yes" if c["agree"] else "**no**"))
    L.append("")
    if rep.get("campaign"):
        c = rep["campaign"]
        if c.get("error"):
            L.append("Campaign rules: **not compared** (%s)." % c["error"])
        else:
            L.append("Campaign rules (`rules_per_model`, folders left out, order of the runs): %d values compared, "
                     "%d disagreements; %s." % (c["n_compared"], c["n_disagree"], "agree" if c["agree"] else "**do not agree**"))
        L.append("")
    worst = rep.get("largest_differences", {})
    if worst:
        L.append("Largest differences found, over all values compared: " +
                 "; ".join("%s %.6g" % (k, v) for k, v in sorted(worst.items())) + ".")
        L.append("")
    any_dis = False
    for r in rep["runs"]:
        c = r.get("comparison")
        if not c:
            continue
        thr = [(n, d) for n, f in c["files"].items() for d in f["threshold_frames"]] + \
              [("summary.json", d) for d in c["summary"]["threshold"]]
        if c["agree"] and not thr:
            continue
        any_dis = any_dis or not c["agree"]
        L.append("## %s" % r["label"])
        L.append("")
        for name, f in c["files"].items():
            for d in f["disagreements"][:12]:
                L.append("- disagreement, `%s`: %s" % (name, json.dumps(d)))
            if f["n_disagree"] > 12:
                L.append("- ... %d disagreements in `%s` in all" % (f["n_disagree"], name))
        for d in c["summary"]["disagreements"][:40]:
            L.append("- disagreement, `summary.json`: %s" % json.dumps(d))
        if c["summary"]["only_primary"]:
            L.append("- keys only in the first program: %s" % ", ".join(c["summary"]["only_primary"][:30]))
        if c["summary"]["only_recheck"]:
            L.append("- keys only in the second program: %s" % ", ".join(c["summary"]["only_recheck"][:30]))
        for name, d in thr[:30]:
            L.append("- threshold frame or what follows from one (listed, both readings; no disagreement), `%s`: %s"
                     % (name, json.dumps(d)))
        if c["summary"]["n_warnings_primary"] != c["summary"]["n_warnings_recheck"]:
            L.append("- warnings of the first program: %s" % json.dumps(c["summary"]["warnings_primary"]))
            L.append("- warnings of the second program: %s" % json.dumps(c["summary"]["warnings_recheck"]))
        L.append("")
    if not any_dis:
        L.append("No disagreement to list.")
        L.append("")
    labels = {}
    for r in rep["runs"]:
        c = r.get("comparison")
        for l in (c["summary"].get("labels_differ", []) if c else []):
            k = l["key"].split("/")
            key = (k[-1] if k[0] == "stats" else l["key"], str(l["primary"]), str(l["recheck"]))
            labels[key] = labels.get(key, 0) + 1
    if labels:
        L.append("Texts that the specification does not fix and that differ (not counted): " +
                 "; ".join("`%s` '%s' against '%s' (%d times)" % (k[0], k[1], k[2], n) for k, n in sorted(labels.items())) + ".")
        L.append("")
    return "\n".join(L) + "\n"


def primary_script():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "analyse_md.py")


def call_primary(arguments, log):
    cmd = [sys.executable, primary_script()] + arguments
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    log.append({"command": " ".join(cmd), "returncode": p.returncode, "tail": p.stdout[-600:]})
    return p.returncode == 0


def spec_version_of(summary):
    """The version of the specification named in the key `spec` of a summary.json, as text, or None."""
    m = re.search(r"version\s+(\d+(?:\.\d+)*)", str((summary or {}).get("spec")))
    return m.group(1) if m else None


def primary_is_stale(run_dir, primary_dir):
    """The first program's output is missing, older than a trajectory file or the manifest of the run,
    or written after another version of the specification than this program follows."""
    s = os.path.join(primary_dir, "summary.json")
    if not os.path.isfile(s):
        return True
    try:
        with open(s) as fh:
            if spec_version_of(json.load(fh)) != SPEC_VERSION:
                return True
    except (ValueError, OSError):
        return True
    t = os.path.getmtime(s)
    return any(os.path.getmtime(os.path.join(run_dir, f)) > t for f in os.listdir(run_dir)
               if f.endswith(".dcd") or f == "manifest.json")


def do_compare(run_dirs, args, top=None):
    """Analyse and compare a list of run folders; with `top` also the campaign rules."""
    log = []
    runs = []
    summaries = {}
    out_root = os.path.abspath(args.out) if getattr(args, "out", None) else None
    for d in run_dirs:
        d = os.path.abspath(d)
        rel = os.path.relpath(d, top) if top else os.path.basename(d)
        a = argparse.Namespace(**vars(args))
        if top:
            a.start = a.model = a.planned_ns = a.system = None
            a.out = os.path.join(out_root, rel, "analysis_recheck") if out_root else None
        else:
            a.out = out_root
        primary_dir = args.primary if (getattr(args, "primary", None) and not top) else os.path.join(d, "analysis")
        entry = {"run_dir": d, "label": rel}
        try:
            summ, odir = analyse_run(d, a, quiet=True)
        except Stop as err:
            entry["error"] = "the second program stopped: %s" % err
            runs.append(entry)
            continue
        summaries[d] = (summ, odir)
        if getattr(args, "run_primary", False) and primary_is_stale(d, primary_dir):
            extra = []
            for opt, val in (("--system", args.system), ("--start", args.start), ("--model", args.model),
                             ("--planned-ns", args.planned_ns)):
                if val is not None and not top:
                    extra += [opt, str(val)]
            call_primary([d] + extra, log)
        entry.update({"n_frames_solute": summ["run"]["n_frames_solute"], "n_frames_full": summ["run"]["n_frames_full"],
                      "role": summ["run"]["role"], "verdict": summ["orientation"]["verdict"]})
        c = compare_run(primary_dir, odir)
        if c.get("error"):
            entry["error"] = c["error"]
        else:
            entry["comparison"] = c
        runs.append(entry)
        print("%-28s %s" % (rel, entry.get("error") or ("agree" if c["agree"] else "DISAGREE") +
                            " (%d values, %d disagreements, %d threshold items, %d keys in one program only)"
                            % (c["n_values"], c["n_disagree"], c["n_threshold"], c["n_keys_one_side"])))
    rep = {"what": top or ", ".join(run_dirs), "spec": "md/analysis/SPEC.md version " + SPEC_VERSION,
           "script_sha256": sha256_file(os.path.abspath(__file__)),
           "primary_script_sha256": sha256_file(primary_script()) if os.path.isfile(primary_script()) else None,
           "runs": runs, "primary_calls": log}
    camp = None
    if top:
        a = argparse.Namespace(**vars(args))
        cdir = os.path.join(out_root, "analysis_campaign_recheck") if out_root else None
        rec, cdir = run_campaign(top, a, summaries=summaries, out_dir=cdir)
        pj = os.path.join(top, "analysis_campaign", "campaign.json")
        if getattr(args, "run_primary", False) and (not os.path.isfile(pj) or any(
                os.path.getmtime(os.path.join(d, "analysis", "summary.json")) > os.path.getmtime(pj)
                for d in summaries if os.path.isfile(os.path.join(d, "analysis", "summary.json")))):
            call_primary(["--campaign", top, "--reuse"], log)
        camp = compare_campaign(pj, rec)
        rep["campaign"] = camp
        rep["campaign_recheck"] = {"rules_per_model": rec["rules_per_model"], "skipped": rec["skipped"],
                                   "failed": rec["failed"], "folder": cdir}
        report_dir = cdir
    else:
        report_dir = out_root or os.path.join(os.path.abspath(run_dirs[0]), "analysis_recheck")
    worst = {}
    for r in runs:
        c = r.get("comparison")
        if not c:
            continue
        for k, v in c["summary"]["max_diff"].items():
            worst["summary " + k] = max(worst.get("summary " + k, 0.0), v)
        for f in c["files"].values():
            for col, st in f["columns"].items():
                if "max_abs_diff" in st:
                    worst["series " + st["kind"]] = max(worst.get("series " + st["kind"], 0.0), st["max_abs_diff"])
    rep["largest_differences"] = worst
    n_err = sum(1 for r in runs if r.get("error"))
    n_dis = sum(1 for r in runs if r.get("comparison") and not r["comparison"]["agree"])
    camp_bad = bool(camp and camp.get("agree") is False)
    camp_missing = bool(camp and camp.get("agree") is None)
    rep["n_runs"] = len(runs)
    rep["n_runs_agree"] = sum(1 for r in runs if r.get("comparison") and r["comparison"]["agree"])
    rep["n_runs_disagree"] = n_dis
    rep["n_runs_not_compared"] = n_err
    rep["agree"] = (n_dis == 0 and n_err == 0 and not camp_bad and not camp_missing)
    # a verdict, a part of a gate, a rule or a role that differs between the two programs through threshold
    # frames or a mean at a limit: allowed by section 14, reported with both readings, and to be read
    to_read = [{"run": r["label"], "keys": r["comparison"]["keys_that_depend_on_threshold_frames"]}
               for r in runs if r.get("comparison") and r["comparison"]["keys_that_depend_on_threshold_frames"]]
    rep["verdicts_that_differ_through_threshold_frames"] = to_read
    rep["n_threshold_items"] = sum(r["comparison"]["n_threshold"] for r in runs if r.get("comparison"))
    if rep["agree"] and to_read:
        rep["result_text"] = ("the two programs agree within section 14 on all %d run(s)%s, BUT in %d run(s) a verdict, "
                              "a gate part, a rule or a role differs through frames at a threshold; both readings "
                              "are listed and must be read"
                              % (len(runs), " and on the campaign rules" if camp else "", len(to_read)))
    elif rep["agree"]:
        rep["result_text"] = "the two programs agree on all %d run(s)%s" % (len(runs), " and on the campaign rules" if camp else "")
    else:
        rep["result_text"] = ("%d of %d run(s) agree, %d disagree, %d could not be compared%s"
                              % (rep["n_runs_agree"], len(runs), n_dis, n_err,
                                 ("; campaign rules: " + ("do not agree" if camp_bad else "not compared" if camp_missing else "agree")) if camp else ""))
    os.makedirs(report_dir, exist_ok=True)
    write_json(os.path.join(report_dir, "recheck_report.json"), rep)
    with open(os.path.join(report_dir, "recheck_report.md"), "w") as fh:
        fh.write(report_markdown(rep))
    print("Result: %s." % rep["result_text"])
    print("Report: %s" % os.path.join(report_dir, "recheck_report.md"))
    rep_status = 1 if (n_dis or camp_bad) else 2 if (n_err or camp_missing) else 3 if to_read else 0
    print("Exit status %d (0 agree; 1 disagreement beyond tolerance; 2 something could not be compared; "
          "3 a verdict differs through threshold frames)." % rep_status)
    return rep_status


# ----------------------------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description="Second, independent implementation of md/analysis/SPEC.md "
                                             "(version %s) and its comparison with analyse_md.py." % SPEC_VERSION)
    ap.add_argument("run_dir", nargs="?", help="run folder to analyse")
    ap.add_argument("--campaign", metavar="DIR", help="analyse the production runs below DIR and apply the rules")
    ap.add_argument("--compare", metavar="RUN_DIR", help="analyse RUN_DIR and compare with analyse_md.py's output")
    ap.add_argument("--compare-campaign", metavar="DIR",
                    help="analyse and compare every run below DIR and the campaign rules (the one command of stage S6)")
    ap.add_argument("--primary", metavar="DIR", help="with --compare: folder of analyse_md.py's output "
                                                     "(default RUN_DIR/analysis)")
    ap.add_argument("--run-primary", action="store_true",
                    help="with --compare / --compare-campaign: call analyse_md.py where its output is missing or old")
    ap.add_argument("--system", metavar="DIR")
    ap.add_argument("--start", choices=["X", "N", "R"])
    ap.add_argument("--model")
    ap.add_argument("--planned-ns", type=float)
    ap.add_argument("--solute-interval-ps", type=float)
    ap.add_argument("--full-interval-ps", type=float)
    ap.add_argument("--out", metavar="DIR", help="output folder (default: RUN_DIR/analysis_recheck, "
                                                 "DIR/analysis_campaign_recheck)")
    ap.add_argument("--no-hash-trajectories", action="store_true")
    ap.add_argument("--skip-name-check", action="store_true")
    args = ap.parse_args(argv)
    modes = [bool(args.run_dir), bool(args.campaign), bool(args.compare), bool(args.compare_campaign)]
    if sum(modes) != 1:
        ap.error("give exactly one of RUN_DIR, --campaign DIR, --compare RUN_DIR, --compare-campaign DIR")
    try:
        if args.run_dir:
            analyse_run(args.run_dir, args)
            return 0
        if args.campaign:
            rec, odir = run_campaign(args.campaign, args)
            for m, rl in rec["rules_per_model"].items():
                print("model %s: %d runs, outcome %s%s" % (m, rl["n_runs"], rl["outcome"], "" if rl["final"] else " (not final)"))
            print("skipped: %s" % ", ".join(s["folder"] for s in rec["skipped"]) if rec["skipped"] else "skipped: none")
            for f in rec["failed"]:
                print("FAILED %s: %s" % (f["folder"], f["error"]))
            print("written to %s" % odir)
            return 0
        if args.compare:
            return do_compare([args.compare], args)
        top = os.path.abspath(args.compare_campaign)
        dirs = find_run_folders(top)
        if not dirs:
            raise Stop("no run folder with a solute_segNNNN.dcd below %s" % top)
        return do_compare(dirs, args, top=top)
    except Stop as err:
        print("STOPPED: %s" % err, file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
