#!/usr/bin/env python3
"""Writes md/gpu/tasks_production.json and md/gpu/tasks_sensitivity.json for a given primary zinc set.
usage: make_tasks.py PRIMARY(2e|2s) OUTDIR  "TEXT ON THE OTHER SET" """
import json, sys
from pathlib import Path
prim, outdir, other_text = sys.argv[1], Path(sys.argv[2]), sys.argv[3]
other = {"2e": "2s", "2s": "2e"}[prim]
name = {"2e": "empirical (MCPB.py step 2e)", "2s": "Seminario (MCPB.py step 2s)"}
env = {"MD": "/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate/md",
       "PY": "/root/cqc-sfa-tools/md/bin/python",
       "SYS": "/root/cqc-sfa-md/systems",
       "RUNS": "/root/cqc-sfa-md/runs"}
LIMIT = 21600


def run_cmd(model, start, rid):
    return (f'"$PY" "$MD/run_md.py" --system "$SYS/{model}/{start}" --out "$RUNS/{model}/{rid}" --run-id {rid} '
            f'--platform CUDA --stage all --prod-ns 100 --segment-ns 5 --resume')


def test(model, rid):
    return (f'python3 "$MD/gpu/queue.py" json "$RUNS/{model}/{rid}/manifest.json" run_id={rid} status=complete '
            f'production.complete=true production.completed_ns=100 test_only=false')


# ---- the nine production runs
M = f"bonded_{prim}"
tasks, before = [], [f"X0_{prim}", "preflight_speed", "preflight_interrupt"]
for rep in (1, 2, 3):
    for start in ("X", "N", "R"):
        rid = f"{start}{rep}"
        tasks.append({"id": rid,
                      "title": f"production run {rid}: 100 ns from start {start}, replica {rep}, zinc set {prim} (primary), terms at N2 as built",
                      "cwd": "$RUNS", "resources": "gpu", "depends_on": before, "stop_grace_s": 120, "timeout_s": LIMIT,
                      "command": run_cmd(M, start, rid), "success_test": test(M, rid)})
        before = [rid]
prod = {"comment": (f"The nine production runs (md/PLAN.md stage S5; entry of 6 Oct 2026, 9:15 PM IST, point 2; md/PRODUCTION.md) "
                    f"with the primary zinc set {prim}, {name[prim]}, systems $SYS/{M}/X, N, R, terms at N2 as built. "
                    f"One after another in the order X1 N1 R1, X2 N2 R2, X3 N3 R3: each waits for the one before it, so the chain "
                    f"stops at the first run that fails; X1 waits for the validation run X0_{prim}, whose gate G-V was judged on "
                    f"the Pi before this file was added. Each run: full protocol of md/protocol.json, 100 ns in segments of 5 ns, "
                    f"--resume, so a run that was stopped goes on from its last complete segment. Time limit of one attempt: "
                    f"{LIMIT // 3600} hours (a run takes about 3 hours); a run that reaches it is ended and continued by the next attempt."),
        "env": env, "tasks": tasks}

# ---- the two sensitivity runs of the rotated start
MO, MN = f"bonded_{other}", f"bonded_{prim}_n2fit"
t_zn = (f"SENSITIVITY RUN of the rotated start, not one of the nine runs and in no count of rule 3 (md/PLAN.md, entry of 6 Oct 2026, 9:15 PM IST, point 3, "
        f"md/PRODUCTION.md). Zinc set {other}, {name[other]}, in place of the primary set {prim}. Terms at N2 as built. {other_text} "
        f"Start R, 100 ns, run id R1, so the same seeds as run R1 of the primary set.")
t_n2 = (f"SENSITIVITY RUN of the rotated start, not one of the nine runs and in no count of rule 3 (md/PLAN.md, entry of 6 Oct 2026, 9:15 PM IST, point 3, "
        f"md/PRODUCTION.md). Zinc set {prim} as in the nine runs. Terms at N2 from md/qm/n2/n2_terms_new_charges_as_built.frcmod, "
        f"the candidate that gate G-R2 did NOT ACCEPT (md/N2.md). It is used only to bracket the terms as built. "
        f"Start R, 100 ns, run id R1, so the same seeds as run R1 of the primary set.")
for t in (t_zn, t_n2, other_text):
    assert not any(c in t for c in '"$`\\!'), t


def marked(model, text):
    d = f'"$RUNS/{model}/R1"'
    return f'mkdir -p {d} && echo "{text}" > "$RUNS/{model}/R1/SENSITIVITY_RUN.txt" && ' + run_cmd(model, "R", "R1")


sens = {"comment": (f"The two sensitivity runs of the rotated start (md/PLAN.md, rule 4 and entry of 6 Oct 2026, 9:15 PM IST, point 3; "
                    f"md/PRODUCTION.md), each once, 100 ns, after the nine production runs: R_zn waits for R3, R_n2 for R_zn. "
                    f"R_zn: the other zinc set ({other}), system $SYS/{MO}/R, folder $RUNS/{MO}/R1. "
                    f"R_n2: the primary zinc set ({prim}) with the candidate terms at N2 that gate G-R2 did not accept, system "
                    f"$SYS/{MN}/R, folder $RUNS/{MN}/R1. Neither enters a count of rule 3. Each task writes the file "
                    f"SENSITIVITY_RUN.txt into its run folder before the run starts. Run id R1 in both, so the seeds are those of run "
                    f"R1 of the primary set. Time limit of one attempt: {LIMIT // 3600} hours."),
        "env": env,
        "tasks": [
            {"id": "R_zn", "title": f"sensitivity run (zinc set {other}): 100 ns from start R, terms at N2 as built; no count of rule 3",
             "cwd": "$RUNS", "resources": "gpu", "depends_on": ["R3"], "stop_grace_s": 120, "timeout_s": LIMIT,
             "command": marked(MO, t_zn), "success_test": test(MO, "R1")},
            {"id": "R_n2", "title": f"sensitivity run (terms at N2: candidate not accepted by gate G-R2; zinc set {prim}): 100 ns from start R; no count of rule 3",
             "cwd": "$RUNS", "resources": "gpu", "depends_on": ["R_zn"], "stop_grace_s": 120, "timeout_s": LIMIT,
             "command": marked(MN, t_n2), "success_test": test(MN, "R1")}]}
for fn, doc in (("tasks_production.json", prod), ("tasks_sensitivity.json", sens)):
    with open(outdir / fn, "w") as f:
        json.dump(doc, f, indent=1)
        f.write("\n")
    print("written", outdir / fn, len(doc["tasks"]), "tasks")
