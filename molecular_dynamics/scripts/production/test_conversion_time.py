#!/usr/bin/env python3
"""Test of md/production/conversion_time.py on made-up series with a known answer.

Each case writes a run folder with analysis/ and analysis_recheck/ (the two CSV files with a few columns in the
order of the real files, and a summary.json with the keys the script reads) into a temporary folder, runs the
script's own functions on it and compares with the answer worked out by hand in the comment of the case.

    python md/production/test_conversion_time.py          prints one line per case; exit status 0 if all pass
"""
import json
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import conversion_time as ct  # noqa: E402

EQ_HEADER = ["frame", "stage", "extension", "time_ns", "lig_rmsd_A", "d1_N1_F131CZ_A", "d2_N1_N67OD1_A", "state",
             "d_N2_T199OG1_A"]
PR_HEADER = ["frame", "segment", "time_ns", "lig_rmsd_A", "d1_N1_F131CZ_A", "d2_N1_N67OD1_A", "state",
             "state_assigned", "d_N2_T199OG1_A"]


def expand(spec):
    """'R*6 I*1 C*123' -> list of states."""
    out = []
    for part in spec.split():
        s, n = part.split("*")
        out += [s] * int(n)
    return out


def write_program(folder, start, start_state, eq_steps, prod_states):
    """eq_steps: list of (stage, extension, length in ps, states of its frames)."""
    os.makedirs(folder)
    total_ps = sum(length for _, _, length, _ in eq_steps)
    rows = []
    frame = 0
    s0 = 0.0
    for stage, ext, length, states in eq_steps:
        assert len(states) == int(length // 10), (stage, ext, len(states), length)
        for k, st in enumerate(states):
            t = (s0 + (k + 1) * 10.0 - total_ps) / 1000.0
            rows.append([frame, stage, ext, "%.6f" % t, "1.0000", "5.0000", "5.0000", st, "2.9000"])
            frame += 1
        s0 += length
    if eq_steps:
        with open(os.path.join(folder, ct.EQUIL_FILE), "w") as fh:
            fh.write(",".join(EQ_HEADER) + "\n")
            for r in rows:
                fh.write(",".join(str(x) for x in r) + "\n")
    with open(os.path.join(folder, ct.PROD_FILE), "w") as fh:
        fh.write(",".join(PR_HEADER) + "\n")
        for n, st in enumerate(prod_states):
            fh.write(",".join(str(x) for x in [n, n // 500 + 1, "%.3f" % ((n + 1) * 0.01), "1.0000", "5.0000",
                                                "5.0000", st, "X", "2.9000"]) + "\n")
    summ = {"run": {"start": start, "model": "made_up", "role": "production", "solute_interval_ps": 10.0},
            "orientation": {"start_state": start_state, "opposite_state": ct.OPPOSITE[start_state],
                            "verdict": "NOT COMPUTED IN THE TEST", "first_passage_ns": None,
                            "first_passage_frame": None},
            "equilibration": {"total_ps": total_ps,
                              "stages": [{"tag": "eq%02d%s" % (s, ("x%d" % e) if e else "")} for s, e, _, _ in eq_steps]}}
    with open(os.path.join(folder, "summary.json"), "w") as fh:
        json.dump(summ, fh)


def held(state="R", extension_steps=0, ext_state=None):
    """Stages 1 to 3 (100 + 300 + 300 ps) and extension steps of 500 ps, all before the release."""
    steps = [(1, 0, 100.0, [state] * 10), (2, 0, 300.0, [state] * 30), (3, 0, 300.0, [state] * 30)]
    for m in range(1, extension_steps + 1):
        steps.append((3, m, 500.0, [ext_state or state] * 50))
    return steps


def released(spec):
    """Stages 4, 5, 6 (400 + 400 + 500 ps = 130 frames) from one string of 130 states."""
    st = expand(spec)
    assert len(st) == 130, len(st)
    return [(4, 0, 400.0, st[:40]), (5, 0, 400.0, st[40:80]), (6, 0, 500.0, st[80:])]


CASES = []


def case(name, start, start_state, eq_first, prod_first, expect, eq_second=None, prod_second=None,
         expect_same=True, comment="", released_frames=130, expect_words=None):
    CASES.append(dict(name=name, start=start, start_state=start_state, eq=(eq_first, eq_second or eq_first),
                      prod=(expand(prod_first), expand(prod_second or prod_first)), expect=expect,
                      expect_same=expect_same, comment=comment, released_frames=released_frames,
                      expect_words=expect_words))


# 1. Start R. Released frames: R in joined frames 0 to 5, one intermediate frame (6), C from frame 7 on. The
#    intermediate frame belongs to R (the state last visited), so the C visit begins at joined frame 7, that is
#    80 ps after the release, in the equilibration (stage 4, frame 70 + 7 = 77 of the equilibration file).
case("converts 80 ps after release, in the equilibration", "R", "R",
     held("R") + released("R*6 I*1 C*123"), "C*1000",
     dict(in_equilibration=True, joined_index=7, frame_in_its_file=77, stage=4, time_from_release_ps=80.0,
          time_in_production_ns=None))

# 2. Start R. The whole released equilibration is R. Production: R in frames 0 to 249, a C visit of 99 frames
#    (250 to 348; one frame short of 1 ns), R in 349 to 358, C from 359 on (641 frames). The first long C visit
#    begins at production frame 359: 3.60 ns in the production, 3.60 + 1.30 = 4.90 ns from the release.
case("converts in the production after a short visit", "R", "R",
     held("R") + released("R*130"), "R*250 C*99 R*10 C*641",
     dict(in_equilibration=False, joined_index=130 + 359, frame_in_its_file=359, stage=None,
          time_from_release_ps=4900.0, time_in_production_ns=3.60))

# 3. Start R. A C visit that begins in stage 6, 30 frames before time zero (joined frames 100 to 129), and goes
#    on for 80 production frames, then R to the end. Counted on the production alone it has 80 frames and is
#    short; on the joined series it has 110 frames and is long. Conversion at joined frame 100, 1010 ps after
#    the release, in the equilibration (stage 6, frame 70 + 100 = 170 of the equilibration file).
case("a visit across time zero is long only on the joined series", "R", "R",
     held("R") + released("R*100 C*30"), "C*80 R*920",
     dict(in_equilibration=True, joined_index=100, frame_in_its_file=170, stage=6, time_from_release_ps=1010.0,
          time_in_production_ns=None))

# 4. Start X (state C). Intermediate frames after C frames (they belong to C) and two R visits of 99 frames
#    each, each ended by a C frame: no long visit in R, so no time.
#    (The first form of this case had "R*99 I*50" and expected no time. That was a mistake of the hand answer,
#    not of the script: intermediate frames belong to the state last visited, so that stretch is one R visit of
#    149 frames. It is kept as case 11.)
case("no conversion: only short visits and intermediate frames", "X", "C",
     held("C") + released("C*60 I*20 C*50"), "C*100 R*99 C*1 R*99 C*1 I*50 C*650", None)

# 5. Start X, mirrored: converts to R in the production. C in frames 0 to 499, intermediate in 500 to 519 (they
#    belong to C), R from 520 on (480 frames): 5.21 ns in the production, 6.51 ns from the release.
case("crystal start that converts in the production (mirrored)", "X", "C",
     held("C") + released("C*130"), "C*500 I*20 R*480",
     dict(in_equilibration=False, joined_index=130 + 520, frame_in_its_file=520, stage=None,
          time_from_release_ps=6510.0, time_in_production_ns=5.21))

# 6. Start R with one extension step of stage 3 (500 ps, 50 frames). Its frames are written as C here although
#    the ligand is held: they come before the release and must not count. Released frames: R in joined frames 0
#    to 19, C from 20 on: 210 ps after the release; frame 70 + 50 + 20 = 140 of the equilibration file.
case("frames of an extension step of stage 3 do not count", "R", "R",
     held("R", extension_steps=1, ext_state="C") + released("R*20 C*110"), "C*1000",
     dict(in_equilibration=True, joined_index=20, frame_in_its_file=140, stage=4, time_from_release_ps=210.0,
          time_in_production_ns=None))

# 7. Start R, intermediate frames from the very first released frame, then C: the intermediate frames before any
#    C or R frame have the state of the start (R), so the C visit begins at the first C frame, joined frame 5:
#    60 ps after the release.
case("intermediate frames at the release belong to the state of the start", "R", "R",
     held("R") + released("I*5 C*125"), "C*500",
     dict(in_equilibration=True, joined_index=5, frame_in_its_file=75, stage=4, time_from_release_ps=60.0,
          time_in_production_ns=None))

# 8. A C visit of exactly 100 frames is long (at least 100), one of 99 is not. Start R, production: R*10 C*100
#    R*890 -> conversion at production frame 10 (0.11 ns; 1.41 ns from the release).
case("a visit of exactly 100 frames is long", "R", "R",
     held("R") + released("R*130"), "R*10 C*100 R*890",
     dict(in_equilibration=False, joined_index=130 + 10, frame_in_its_file=10, stage=None,
          time_from_release_ps=1410.0, time_in_production_ns=0.11))

# 9. The two programs differ in one state: in the second program the frame before the visit is C as well, so
#    its visit begins one frame earlier. The script must say that the frames differ and report no time.
case("the two programs give different frames", "R", "R",
     held("R") + released("R*130"), "R*300 C*700",
     dict(in_equilibration=False, joined_index=130 + 300, frame_in_its_file=300, stage=None,
          time_from_release_ps=4310.0, time_in_production_ns=3.01),
     prod_second="R*299 C*701", expect_same=False)

# 10. Start R, the run goes to C for 1.5 ns and comes back: the conversion time is that of the first long visit
#     in C whatever happens later (the verdict of the run is another matter). R*200 C*150 R*650 -> frame 200:
#     2.01 ns in the production, 3.31 ns from the release.
case("a long visit that does not last still gives the time", "R", "R",
     held("R") + released("R*130"), "R*200 C*150 R*650",
     dict(in_equilibration=False, joined_index=130 + 200, frame_in_its_file=200, stage=None,
          time_from_release_ps=3310.0, time_in_production_ns=2.01))

# 11. Intermediate frames lengthen the visit they follow (section 4, point 1: "an intermediate frame belongs to
#     the state last visited"). Start X: C in frames 0 to 99, R in 100 to 159 (60 frames), intermediate in 160 to
#     199 (40 frames), then C. The R visit has 60 + 40 = 100 assigned frames and is long: 1.01 ns in the
#     production, 2.31 ns from the release.
case("intermediate frames count to the visit they follow", "X", "C",
     held("C") + released("C*130"), "C*100 R*60 I*40 C*800",
     dict(in_equilibration=False, joined_index=130 + 100, frame_in_its_file=100, stage=None,
          time_from_release_ps=2310.0, time_in_production_ns=1.01))

# 12. A folder without an equilibration series (no run of the campaign has this; audit of 8 Oct 2026). Start R,
#     production R*50 C*200: the first long C visit begins at production frame 50, 0.51 ns in the production.
#     The release is not on record, so there is no time from the release, and the wording must say so instead
#     of the program ending with an error.
case("no equilibration series: the production time alone, with a plain message", "R", "R",
     [], "R*50 C*200",
     dict(in_equilibration=False, joined_index=50, frame_in_its_file=50, stage=None,
          time_from_release_ps=None, time_in_production_ns=0.51),
     released_frames=0,
     expect_words="0.51 ns in the production (no equilibration series in the folder: the time from the release "
                  "is not known)")


def run_tests():
    n_ok = 0
    with tempfile.TemporaryDirectory() as tmp:
        for i, c in enumerate(CASES, 1):
            d = os.path.join(tmp, "made_up", "T%d" % i)
            write_program(os.path.join(d, "analysis"), c["start"], c["start_state"], c["eq"][0], c["prod"][0])
            write_program(os.path.join(d, "analysis_recheck"), c["start"], c["start_state"], c["eq"][1], c["prod"][1])
            rec = ct.analyse_run(d, "made_up/T%d" % i)
            first = rec["programs"]["analysis"]["conversion"]
            problems = []
            if rec["same_frame"] != c["expect_same"]:
                problems.append("same_frame %s, expected %s" % (rec["same_frame"], c["expect_same"]))
            exp = c["expect"]
            if exp is None:
                if first is not None:
                    problems.append("a conversion was found: %s" % first)
                if c["expect_same"] and not rec["reported_as"].startswith("none"):
                    problems.append("wording %r" % rec["reported_as"])
            else:
                if first is None:
                    problems.append("no conversion found")
                else:
                    for k, v in exp.items():
                        got = first[k]
                        same = (abs(got - v) < 1e-6) if isinstance(v, float) and got is not None else (got == v)
                        if not same:
                            problems.append("%s = %r, expected %r" % (k, got, v))
                if c["expect_same"]:
                    if rec["conversion"] is None:
                        problems.append("the run record holds no conversion")
                else:
                    if rec["conversion"] is not None or "NOT REPORTED" not in rec["reported_as"]:
                        problems.append("a time is reported although the programs differ")
            if c["expect_words"] is not None and rec["reported_as"] != c["expect_words"]:
                problems.append("wording %r" % rec["reported_as"])
            # the released part must be 130 frames = 1300 ps in every case with an equilibration series
            for sub in ("analysis", "analysis_recheck"):
                p = rec["programs"][sub]
                if c["released_frames"] == 0:
                    if p["n_released_equilibration_frames"] != 0 or p["released_part_of_equilibration_ps"] is not None:
                        problems.append("%s: a released part without an equilibration series" % sub)
                elif p["n_released_equilibration_frames"] != 130 or abs(p["released_part_of_equilibration_ps"] - 1300.0) > 1e-6:
                    problems.append("%s: released part %s frames, %s ps" % (sub, p["n_released_equilibration_frames"],
                                                                         p["released_part_of_equilibration_ps"]))
            ok = not problems
            n_ok += ok
            print("%-4s case %2d: %-66s %s" % ("ok" if ok else "FAIL", i, c["name"],
                                              rec["reported_as"] if ok else "; ".join(problems)))
        # the function on the plain rule, against a slow reading of section 4 written out again here
        import random
        rnd = random.Random(20261008)
        n_rand = 0
        for trial in range(300):
            start_state = rnd.choice("CR")
            states = []
            while len(states) < 1500:
                states += [rnd.choice("CCRRI")] * rnd.choice([1, 2, 5, 40, 99, 100, 101, 250])
            states = states[:1500]
            got, _, _ = ct.first_long_visit_opposite(states, start_state)
            # slow reading: walk the frames, keep the assigned state, measure each stretch
            a = []
            last = start_state
            for s in states:
                last = s if s in "CR" else last
                a.append(last)
            want = None
            k = 0
            while k < len(a):
                m = k
                while m < len(a) and a[m] == a[k]:
                    m += 1
                if want is None and a[k] != start_state and m - k >= 100:
                    want = k
                k = m
            n_rand += (got == want)
        print("%-4s random series against a second reading of section 4: %d of 300 the same"
              % ("ok" if n_rand == 300 else "FAIL", n_rand))
    total = len(CASES) + 1
    n_ok += (n_rand == 300)
    print("%d of %d checks passed" % (n_ok, total))
    return 0 if n_ok == total else 1


if __name__ == "__main__":
    sys.exit(run_tests())
