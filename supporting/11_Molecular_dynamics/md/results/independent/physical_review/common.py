"""Shared loader for the review scripts (reads the npz files written by extract.py)."""
import os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
RUNLIST = [("bonded_2e", r) for r in ("X1", "X2", "X3", "N1", "N2", "N3", "R1", "R2", "R3")] + [("bonded_2s", "R1"), ("bonded_2e_n2fit", "R1")]
def label(model, run):
    return {"bonded_2e": "2e", "bonded_2s": "2s", "bonded_2e_n2fit": "n2fit"}[model] + "/" + run
class Run:
    def __init__(self, model, run):
        d = np.load(f"{HERE}/data_{model}_{run}.npz", allow_pickle=False)
        self.model, self.run, self.name = model, run, label(model, run)
        self.meta = json.loads(str(d["meta"]))
        self.s = {k[2:]: d[k] for k in d.files if k.startswith("s_")}
        self.f = {k[2:]: d[k] for k in d.files if k.startswith("f_")}
        self.start = run[0]
        self.t = self.s["t_ps"]; self.tf = self.f["t_ps"]
        self.trel = self.meta["t_release_ps"]
        self.prod = self.t > 0; self.prodf = self.tf > 0
        self.free = self.t > self.trel            # from the release of the ligand (eq04) to the end
        self.freef = self.tf > self.trel
        # state of the solute frame at the time of every full frame
        idx = {round(float(t), 3): i for i, t in enumerate(self.t)}
        self.f2s = np.array([idx[round(float(t), 3)] for t in self.tf])
        self.state = self.s["state"]              # 0 C, 1 R, 2 I
        self.statef = self.state[self.f2s]
def hb(d, a, dcut=3.5, acut=135.0):
    return (d <= dcut) & (a >= acut)
def circ_mean(a_deg):
    a = np.radians(np.asarray(a_deg, float))
    return float(np.degrees(np.arctan2(np.sin(a).mean(), np.cos(a).mean())))
def circ_sd(a_deg):
    a = np.radians(np.asarray(a_deg, float))
    R = np.hypot(np.sin(a).mean(), np.cos(a).mean())
    return float(np.degrees(np.sqrt(max(0.0, -2.0 * np.log(max(R, 1e-12))))))
def assigned(state, start_state):
    out = np.empty(len(state), int); cur = start_state
    for i, s in enumerate(state):
        if s in (0, 1): cur = s
        out[i] = cur
    return out
def visits(a):
    """maximal stretches of equal value: list of (value, first, length)"""
    out = []; i = 0; n = len(a)
    while i < n:
        j = i
        while j + 1 < n and a[j + 1] == a[i]: j += 1
        out.append((int(a[i]), i, j - i + 1)); i = j + 1
    return out
def ms(a, fmt="%.3f"):
    a = np.asarray(a, float)
    return (fmt + " +- " + fmt + " (" + fmt + " to " + fmt + ")") % (a.mean(), a.std(), a.min(), a.max())
