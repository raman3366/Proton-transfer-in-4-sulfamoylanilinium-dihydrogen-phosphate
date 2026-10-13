#!/usr/bin/env python3
"""Build the DFT models M1-M7 from the deposited CIF and write their case files.

M1 ion pair: the asymmetric-unit cation plus the anion translated to (x, y, z-1), i.e. the partner of the
strongest cation->anion bond N1-H1A...O5 (2.760 A, 174 deg), which also closes the R2,2(8) ring with
C6-H6...O4.  Fragments are cut from M1 so every model starts from crystal coordinates:
M2 cation, M3 H2PO4-, M4 neutral sulfanilamide (M2 minus H1A), M5 sulfonamidate (M4 minus H2A),
M6 H3PO4 (M3 plus a proton on O5, placed along O5->H1A), M7 urea (RDKit ETKDG; NLO reference only).
Also the def2-TZVP twins sfa_pair_pcm_d3_tz and sulfanilamide_pcm_d3_tz (basis-set sensitivity), copied from
their parent case files with a case-level "basis".

Usage (from the repo root): python structures/build_models.py
Writes molecules/sfa/input/*.xyz and molecules/sfa/cases/*.json.
"""
import json
from pathlib import Path

import gemmi
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
CIF = ROOT / "e-70-0o997-sup1.cif"
INP = ROOT / "molecules/sfa/input"
CAS = ROOT / "molecules/sfa/cases"

CATION = ["C1", "C2", "C3", "C4", "C5", "C6", "N1", "N2", "O1", "O2", "S1",
          "H2", "H3", "H5", "H6", "H1A", "H1B", "H1C", "H2A", "H2B"]
ANION = ["P1", "O3", "O4", "O5", "O6", "H3A", "H4A"]
COMMON = {"ccdc": 1016950, "source_doi": "10.1107/S1600536814017462", "source_cif": "e-70-0o997-sup1.cif",
          "cell": [4.8041, 10.8564, 10.3862, 90.0, 101.067, 90.0], "spin": 0}


def load():
    ss = gemmi.make_small_structure_from_block(gemmi.cif.read(str(CIF)).sole_block())
    sites = {s.label: s for s in ss.sites}

    def cart(lbl, shift=(0, 0, 0)):
        f = sites[lbl].fract
        p = ss.cell.orthogonalize(gemmi.Fractional(f.x + shift[0], f.y + shift[1], f.z + shift[2]))
        return np.array([p.x, p.y, p.z])

    labels = CATION + [f"{a}#1_554" for a in ANION]
    syms = [sites[l].element.name for l in CATION] + [sites[l].element.name for l in ANION]
    xyz = np.array([cart(l) for l in CATION] + [cart(l, (0, 0, -1)) for l in ANION])
    return labels, syms, xyz


def write_xyz(path, syms, xyz, comment):
    with open(path, "w") as f:
        f.write(f"{len(syms)}\n{comment}\n")
        for s, (x, y, z) in zip(syms, xyz):
            f.write(f"{s:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")


def case(name, xyz_name, title, charge, solvent, disp, labels, extra=None, notes=""):
    c = dict(COMMON, title=title, xyz=f"molecules/sfa/input/{xyz_name}.xyz", charge=charge,
             cif_labels=labels, disp=disp, notes=notes)
    if solvent:
        c["solvent"] = solvent
    c.update(extra or {})
    (CAS / f"{name}.json").write_text(json.dumps(c, indent=1) + "\n")


def basis_sensitivity_cases():
    """def2-TZVP twins of the primary PCM models M1 and M4 (plan P2.6, basis-set sensitivity).

    Each is a copy of its parent case with a case-level "basis" (read by dft/common.mol_for_case) and its own
    name, so results land in results/<parent>_tz/ and nothing of the parent changes.  Same crystal input
    geometry; the burst starts the optimisation from the parent's 6-311++G(d,p) optimum (copied in as
    preopt.xyz, which dft/optimize.py then uses instead of its own small-basis pre-optimisation)."""
    for parent in ("sfa_pair_pcm_d3", "sulfanilamide_pcm_d3"):
        c = json.loads((CAS / f"{parent}.json").read_text())
        c["basis"] = "def2-tzvp"
        c["title"] = (c["title"].split(" -- ")[0].replace("B3LYP-D3(BJ)", "B3LYP-D3(BJ)/def2-TZVP")
                      + " -- basis-set sensitivity check")
        c["notes"] = (c.get("notes", "") + f" | def2-TZVP twin of {parent} (basis-set sensitivity); start the "
                      f"optimisation from results/{parent}/opt.xyz.").strip(" |")
        c["basis_parent_case"] = parent
        (CAS / f"{parent}_tz.json").write_text(json.dumps(c, indent=1) + "\n")


def main():
    INP.mkdir(parents=True, exist_ok=True); CAS.mkdir(parents=True, exist_ok=True)
    labels, syms, xyz = load()
    ix = {l.split("#")[0]: i for i, l in enumerate(labels)}
    d = lambda a, b: float(np.linalg.norm(xyz[ix[a]] - xyz[ix[b]]))
    assert abs(d("N1", "O5") - 2.760) < 0.002, d("N1", "O5")

    # M1 ion pair
    write_xyz(INP / "sfa_pair.xyz", syms, xyz, "4-sulfamoylanilinium dihydrogen phosphate ion pair, CIF coordinates "
              "(anion at x,y,z-1; N1-H1A...O5 2.760 A)")
    pair_extra = {
        "chemical_name": "4-sulfamoylanilinium dihydrogen phosphate", "formula": "C6H11N2O6PS",
        "fragments": {"4-sulfamoylanilinium cation": list(range(20)), "dihydrogen phosphate anion": list(range(20, 27))},
        "hbonds": {"N1-H1A...O5": [ix["N1"], ix["H1A"], ix["O5"]], "C6-H6...O4": [ix["C6"], ix["H6"], ix["O4"]]},
        "distances": {"N1...O5 (X-ray 2.760 A)": [ix["N1"], ix["O5"]], "C6...O4 (X-ray 3.288 A)": [ix["C6"], ix["O4"]]},
        "planes": {"benzene_ring": [ix["C1"], ix["C3"], ix["C5"]], "hbond_ring_R22_8": [ix["N1"], ix["O5"], ix["O4"]],
                   "sulfonamide_S1_N2_O1": [ix["S1"], ix["N2"], ix["O1"]]},
        "proton": {"H": ix["H1A"], "donor_N": ix["N1"], "acceptor_O": ix["O5"]},
    }
    note = ("Ion pair assembled from the CIF: cation of the asymmetric unit and the H2PO4- anion translated by (0,0,-1), "
            "joined by N1-H1A...O5 (2.760 A, 174 deg) and C6-H6...O4 (3.288 A), which together close an R2,2(8) ring. "
            "X-ray H positions are the deposited ones (O-H restrained with DFIX 0.82 A).")
    case("sfa_pair_pcm_d3", "sfa_pair", "SAPH ion pair, B3LYP-D3(BJ), IEF-PCM (water) -- primary model", 0, "water", "d3bj",
         labels, pair_extra, note)
    case("sfa_pair_gas_d3", "sfa_pair", "SAPH ion pair, B3LYP-D3(BJ), gas phase", 0, None, "d3bj", labels, pair_extra, note)
    case("sfa_pair_pcm_b3lyp", "sfa_pair", "SAPH ion pair, plain B3LYP, IEF-PCM (water) -- dispersion sensitivity check",
         0, "water", "none", labels, pair_extra, note)

    # fragments (all B3LYP-D3(BJ); PCM primary, gas for comparison)
    cat = [ix[a] for a in CATION]; an = [ix[a] for a in ANION]
    frag = {}
    frag["cation"] = ([labels[i] for i in cat], [syms[i] for i in cat], xyz[cat], +1, "4-sulfamoylanilinium cation")
    frag["h2po4"] = ([labels[i] for i in an], [syms[i] for i in an], xyz[an], -1, "dihydrogen phosphate anion")
    neu = [i for i in cat if labels[i] != "H1A"]
    frag["sulfanilamide"] = ([labels[i] for i in neu], [syms[i] for i in neu], xyz[neu], 0, "sulfanilamide (neutral)")
    sna = [i for i in neu if labels[i] != "H2A"]
    frag["sulfonamidate"] = ([labels[i] for i in sna], [syms[i] for i in sna], xyz[sna], -1,
                             "sulfanilamide N1-deprotonated sulfonamidate anion (4-aminobenzenesulfonamidate)")
    u = xyz[ix["H1A"]] - xyz[ix["O5"]]; hpos = xyz[ix["O5"]] + 0.97 * u / np.linalg.norm(u)
    frag["h3po4"] = ([labels[i] for i in an] + ["H5P"], [syms[i] for i in an] + ["H"], np.vstack([xyz[an], hpos]), 0,
                     "phosphoric acid (H2PO4- + proton on O5)")
    for k, (lab, sy, xx, q, nm) in frag.items():
        write_xyz(INP / f"{k}.xyz", sy, xx, f"{nm}, cut from the SAPH ion pair (CIF coordinates)")
        envs = [("pcm", "water"), ("gas", None)] if k != "sulfonamidate" else [("pcm", "water")]
        for tag, solv in envs:
            case(f"{k}_{tag}_d3", k, f"{nm}, B3LYP-D3(BJ), {'IEF-PCM (water)' if solv else 'gas phase'}", q, solv, "d3bj",
                 lab, {"chemical_name": nm}, "Fragment cut from the SAPH ion pair; see sfa_pair notes.")

    # M7 urea (NLO reference, gas phase)
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles("NC(N)=O"))
    AllChem.EmbedMolecule(m, randomSeed=7); AllChem.MMFFOptimizeMolecule(m)
    conf = m.GetConformer()
    usy = [a.GetSymbol() for a in m.GetAtoms()]
    uxyz = np.array([list(conf.GetAtomPosition(i)) for i in range(m.GetNumAtoms())])
    write_xyz(INP / "urea.xyz", usy, uxyz, "urea, RDKit ETKDG + MMFF (NLO reference)")
    case("urea_gas_d3", "urea", "urea, B3LYP-D3(BJ), gas phase -- NLO reference", 0, None, "d3bj",
         [f"{s}{i+1}" for i, s in enumerate(usy)], {"chemical_name": "urea"}, "Standard NLO reference molecule.")
    basis_sensitivity_cases()
    print("models:", sorted(p.stem for p in INP.glob("*.xyz")))
    print("cases :", sorted(p.stem for p in CAS.glob("*.json")))


if __name__ == "__main__":
    main()
