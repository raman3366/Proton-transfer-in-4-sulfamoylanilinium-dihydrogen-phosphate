# the check of md/RUN.md, section 4: every file of every step against the sha256 the driver recorded
import hashlib, json, os, sys
def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()
for d in sys.argv[1:]:
    if not os.path.isfile(os.path.join(d, "manifest.json")):
        continue
    m = json.load(open(os.path.join(d, "manifest.json")))
    steps = [s for s in [m["minimisation"]] + m["equilibration"]["stages"] + m["production"]["segments"] if s]
    missing, changed, n = [], [], 0
    for s in steps:
        for name, info in s["files"].items():
            p = os.path.join(d, name); n += 1
            if not os.path.isfile(p):
                missing.append(name)
            elif sha256(p) != info["sha256"]:
                changed.append(name)
    h = m.get("hydration_check", {})
    print(m["run_id"], m["status"], m["production"]["completed_ns"], "ns", "TEST" if m["test_only"] else "production protocol",
          f"{n} files", "checksums OK" if not missing + changed else f"MISSING {len(missing)} CHANGED {len(changed)}: " + " ".join(missing + changed),
          "| hydration check:", h.get("final_decision") or h.get("status"), "| extension steps:", h.get("n_extensions", 0))
