"""
Step 7 -- assemble results/<case>/RESULTS.md: methods, tables, figures and draft interpretation paragraphs
for the six requested analyses. Numbers come from the JSON/CSV files written by the earlier stages.
Usage: python dft/report.py <case>
"""
import sys, json, csv, datetime, numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
from common import *
import pyscf, geometric, matplotlib

case = load_case(sys.argv[1]); d = case['dir']
J = lambda f: json.load(open(d / f)) if (d / f).exists() else None
opt, frq, fmo, mep, pc = J('opt_summary.json'), J('freq.json'), J('fmo.json'), J('mep_surface.json'), J('proton_check.json')
syms, xyz = read_xyz(d / 'opt.xyz'); lab = case.get('cif_labels') or atom_labels(syms)
solv = case.get('solvent')
transferred = bool(pc and pc.get('proton_transferred'))
# after a proton transfer the fragment names in the case file no longer describe the species;
# say so once rather than silently printing the wrong charge state
nm = lambda k: (k + ' (protonation changed)') if transferred and ('cation' in k or 'anion' in k) else k
prot = case.get('proton') or {}
pD = lab[prot['donor_N']] if prot else 'the donor'
pH = lab[prot['H']] if prot else 'H'
pA = lab[prot['acceptor_O']] if prot else 'the acceptor'
frag_list = list(case.get('fragments') or {})
species = ' + '.join(frag_list) if frag_list else case.get('chemical_name', case['name'])
mwfn_ver = ''
for f in (d / 'logs').glob('mwfn_*.out'):
    import re; m = re.search(r'Multiwfn\s+([\d.]+(?:\s*\(dev\))?)', f.read_text(errors='replace')); mwfn_ver = m.group(1) if m else ''; break
fig = lambda name, cap: f'![{cap}](fig/{name})\n\n*{cap}*\n\n' if (d / 'fig' / name).exists() else f'*(figure fig/{name} not generated)*\n\n'
L = []
L.append(f"# DFT results: {case['title']}\n\nGenerated {datetime.datetime.now():%Y-%m-%d %H:%M} by `dft/report.py`. All numbers below are read from the files in this folder; "
         f"the interpretation paragraphs are drafts in the style of comparable papers and should be edited by the chemist.\n")
if case.get('chemical_name'):
    L.append(f"**Compound.** {case['chemical_name']} ({case.get('formula','')}), CCDC {case.get('ccdc','?')}"
             + (f", original structure paper doi:{case['source_doi']}" if case.get('source_doi') else '') + ". "
             + f"Model system: {species} ({len(syms)} atoms).\n")
if case.get('notes'):
    L.append(f"\n**How this model was built.** {case['notes']}\n")
L.append("## Computational methods (draft paragraph)\n")
L.append(f"All calculations were performed with density functional theory at the B3LYP/{basis_description(case)} level using PySCF {pyscf.__version__} "
         f"(B3LYP with the VWN-RPA local correlation, i.e. the same definition as in Gaussian; density fitting with the default auxiliary basis; "
         f"integration grid level 3). Geometry optimisation used the geomeTRIC {geometric.__version__} optimiser with Gaussian's default convergence "
         f"criteria (maximum force 4.5e-4 Ha/Bohr, RMS force 3.0e-4 Ha/Bohr){', in an IEF-PCM continuum of ' + solv if solv else ', in the gas phase'}. "
         f"Harmonic frequencies were computed {('by ' + frq['hessian_method']) if frq else 'analytically'} at the same level to confirm a true minimum. "
         f"Frontier molecular orbitals, the molecular electrostatic potential on the 0.001 a.u. electron-density isosurface, the reduced density gradient "
         f"(NCI), the electron localisation function (ELF) and the localised orbital locator (LOL) were analysed with Multiwfn {mwfn_ver}. "
         f"Isosurfaces were rendered with ParaView 6.1.1 and 2-D maps with matplotlib {matplotlib.__version__}.\n")
# ---- 1. geometry
L.append("## 1. Optimised geometry and frequency check\n")
if opt:
    L.append(f"Final energy: {opt['E_final_Ha']:.6f} Ha after {opt['steps']} optimisation steps ({opt['wall_min']:.0f} min). "
             f"Bond lengths deviate from the X-ray structure by RMSD {opt['bond_rmsd_A']:.3f} Å (max {opt['bond_maxdev_A']:.3f} Å; X-ray H positions are systematically short, so the largest deviations are X-H bonds).\n")
if frq:
    fr = frq['frequencies_cm-1']; L.append(f"Frequencies: {len(fr)} modes, **{frq['n_imaginary']} imaginary** → {'the structure is a true minimum' if frq['n_imaginary']==0 else 'NOT a minimum'}; "
                                            f"lowest modes {', '.join(f'{v:.0f}' for v in fr[:5])} cm⁻¹; ZPE {frq['ZPE_Ha']*HARTREE_EV:.3f} eV.\n")
if pc:
    L.append(f"Proton check after the optimisation ({'IEF-PCM ' + solv if solv else 'gas phase'}): {pH}–{pD} {pc.get('d_H_N_A', float('nan')):.3f} Å, "
             f"{pH}···{pA} {pc.get('d_H_O_A', float('nan')):.3f} Å → "
             f"**{f'the proton migrated from {pD} to {pA}, so the optimised species is the neutral hydrogen-bonded complex rather than the salt' if pc.get('proton_transferred') else 'the proton stays on ' + pD + ', i.e. the ionic (salt) form is retained'}**.\n")
    if pc.get('bonds_broken') or pc.get('bonds_formed'):
        L.append(f"Connectivity change during the optimisation: bonds broken {pc.get('bonds_broken') or 'none'}, formed {pc.get('bonds_formed') or 'none'}.\n")
L.append(fig('structure_labels.png', 'Optimised structure with the atom numbering of the CIF'))
if opt and opt.get('hbonds'):
    L.append("| Hydrogen bond | D–H (Å) | H···A (Å) | D···A (Å) | D–H···A (°) |\n|---|---|---|---|---|\n")
    for k, v in opt['hbonds'].items():
        for tag, g in (('X-ray', v['start']), ('DFT', v['dft'])):
            L.append(f"| {k} ({tag}) | {g['D-H_A']:.3f} | {g['H..A_A']:.3f} | {g['D..A_A']:.3f} | {g['DHA_deg']:.1f} |\n")
    L.append('\n')
s0, x0 = read_xyz(ROOT / case['xyz'])
with open(d / 'geometry_vs_start.csv', 'w') as f:      # (re)write with CIF labels
    f.write('bond,start_A,dft_A,diff_A\n')
    for i, j, d0 in bonds_from_geometry(s0, x0):
        d1 = float(np.linalg.norm(xyz[i] - xyz[j])); f.write(f'{lab[i]}-{lab[j]},{d0:.4f},{d1:.4f},{d1-d0:+.4f}\n')
rows = list(csv.DictReader(open(d / 'geometry_vs_start.csv')))
L.append("Selected bond lengths (Å), heavy atoms only:\n\n| Bond | X-ray | DFT | Δ |\n|---|---|---|---|\n")
for r in rows:
    a, b = r['bond'].split('-')
    if not (a.startswith('H') or b.startswith('H')): L.append(f"| {r['bond']} | {float(r['start_A']):.3f} | {float(r['dft_A']):.3f} | {float(r['diff_A']):+.3f} |\n")
L.append('\nFull list including X–H bonds: `geometry_vs_start.csv`.\n')
# ---- 2. FMO
L.append("## 2. Frontier molecular orbitals and global reactivity descriptors\n")
if fmo:
    L.append(open(d / 'fmo_table.md').read() + '\n')
    fc = fmo.get('fragment_composition')
    if fc:
        L.append("Fragment composition of the frontier orbitals (Mulliken partition, %):\n\n| Orbital | " + ' | '.join(nm(k) for k in next(iter(fc.values())).keys()) + " |\n|---|" + '---|' * len(next(iter(fc.values()))) + "\n")
        for o, comp in fc.items(): L.append(f"| {o} | " + ' | '.join(f'{v:.0f}' for v in comp.values()) + " |\n")
        L.append('\n')
    L.append(fig('mo_levels.png', 'Frontier orbital energy levels') + fig('homo.png', 'HOMO (isovalue ±0.02 a.u.)') + fig('lumo.png', 'LUMO (isovalue ±0.02 a.u.)'))
    hl = ''
    if fc:
        H = max(fc['HOMO'], key=fc['HOMO'].get); Lu = max(fc['LUMO'], key=fc['LUMO'].get)
        hl = (f" The HOMO is localised mainly on the {nm(H)} ({fc['HOMO'][H]:.0f} %) whereas the LUMO resides on the {nm(Lu)} ({fc['LUMO'][Lu]:.0f} %), "
              + ("so the HOMO→LUMO excitation transfers charge between the two partners." if H != Lu else "so the HOMO→LUMO excitation is localised on that fragment."))
    L.append(f"*Draft interpretation.* The HOMO and LUMO energies are {fmo['HOMO_eV']:.3f} and {fmo['LUMO_eV']:.3f} eV, giving an energy gap of {fmo['gap_eV']:.3f} eV.{hl} "
             f"From these, the ionisation potential ({fmo['IP_eV']:.3f} eV) and electron affinity ({fmo['EA_eV']:.3f} eV) yield an electronegativity of {fmo['electronegativity_chi_eV']:.3f} eV, "
             f"a chemical hardness of {fmo['hardness_eta_eV']:.3f} eV and softness of {fmo['softness_S_eV-1']:.3f} eV⁻¹; the electrophilicity index is {fmo['electrophilicity_omega_eV']:.3f} eV. "
             f"The {'moderate' if fmo['gap_eV'] > 3 else 'small'} gap {'indicates kinetic stability combined with' if fmo['gap_eV'] > 3 else 'reflects the ion-pair (charge-transfer) character of the frontier orbitals and suggests'} "
             f"significant polarisability, consistent with the interest in this salt for optical applications. The dipole moment is {fmo['dipole_Debye']:.2f} D.\n")
# ---- 3. MEP
L.append("## 3. Molecular electrostatic potential (MEP)\n")
if mep:
    L.append(f"ESP on the 0.001 a.u. density surface: global minimum {mep['global_min_au']*627.509:.1f} kcal/mol ({mep['global_min_au']:.4f} a.u.), global maximum {mep['global_max_au']*627.509:.1f} kcal/mol ({mep['global_max_au']:.4f} a.u.); "
             f"surface area {mep['area_A2']:.1f} Å² (positive {mep['positive_area_A2']:.1f}, negative {mep['negative_area_A2']:.1f}); molecular polarity index {mep['MPI_eV']:.2f} eV; polar surface fraction {mep['polar_area_pct']:.0f} %.\n")
    ext = list(csv.DictReader(open(d / 'mep_extrema.csv')))
    mins = [e for e in ext if e['type'] == 'minimum'][:4]; maxs = sorted([e for e in ext if e['type'] == 'maximum'], key=lambda e: -float(e['ESP_kcal_mol']))[:4]
    L.append("| Extremum | ESP (kcal/mol) | ESP (a.u.) | Nearest atoms |\n|---|---|---|---|\n")
    for e in mins + maxs: L.append(f"| {e['type']}{' (global)' if e['global'] else ''} | {float(e['ESP_kcal_mol']):.1f} | {float(e['ESP_au']):.4f} | {e['nearest_atoms']} |\n")
    L.append('\n' + fig('mep.png', 'MEP mapped on the 0.001 a.u. electron-density isosurface (red: negative, blue: positive)'))
    minsite = mins[0]['nearest_atoms'].split(',')[0] if mins else '…'
    maxsite = maxs[0]['nearest_atoms'].split(',')[0] if maxs else '…'
    lone_pair_kind = {'O': 'oxygen lone pairs', 'N': 'nitrogen lone pairs', 'Cl': 'chloride lone pairs',
                      'S': 'sulfur lone pairs'}.get(minsite.strip()[0], 'lone pairs')
    L.append(f"*Draft interpretation.* The most negative potential ({mep['global_min_au']*627.509:.1f} kcal/mol) lies near {minsite}, "
             f"identifying the {lone_pair_kind} there as the preferred sites for electrophilic attack and as hydrogen-bond acceptors. "
             f"The most positive potential ({mep['global_max_au']*627.509:.1f} kcal/mol) surrounds {maxsite}, which is therefore the "
             f"nucleophilic-attack and hydrogen-bond-donor site. "
             + (f"After the proton transfer the surface is much less polarised than in the salt form. "
                if transferred else
                f"The strongly polarised surface (polar fraction {mep['polar_area_pct']:.0f} %) reflects the ionic character of the assembly "
                f"and the electrostatic complementarity that drives the hydrogen bonding observed in the crystal. ")
             + f"The molecular polarity index is {mep['MPI_eV']:.2f} eV.\n")
# ---- 4. NCI
L.append("## 4. Reduced density gradient (RDG) / non-covalent interactions\n")
if (d / 'nci_spikes.csv').exists():
    sp = list(csv.DictReader(open(d / 'nci_spikes.csv')))
    L.append("Spikes in the RDG vs sign(λ₂)ρ plot (minimum RDG per bin of sign(λ₂)ρ):\n\n| sign(λ₂)ρ (a.u.) | min RDG | interaction type |\n|---|---|---|\n")
    for r in sp: L.append(f"| {float(r['sign_l2_rho_au']):+.3f} | {float(r['min_RDG']):.2f} | {r['type']} |\n")
    L.append('\n')
L.append(fig('rdg_scatter.png', 'RDG versus sign(λ₂)ρ') + fig('rdg_scatter_wide.png', 'RDG versus sign(λ₂)ρ, extended range showing the strong hydrogen-bond spike') + fig('nci.png', 'NCI isosurface (RDG = 0.5 a.u.) coloured by sign(λ₂)ρ from −0.035 (blue) to +0.02 a.u. (red)'))
hb = [r for r in (list(csv.DictReader(open(d / 'nci_spikes.csv'))) if (d / 'nci_spikes.csv').exists() else []) if 'H-bond' in r['type']]
hbtxt = (f" The hydrogen-bond spike appears at sign(λ₂)ρ = {float(hb[0]['sign_l2_rho_au']):+.3f} a.u. (minimum RDG {float(hb[0]['min_RDG']):.2f}); a value below about −0.05 a.u. indicates an unusually strong, partly covalent hydrogen bond, consistent with the short H···A distance." if hb else '')
if not transferred: L.append("Note: Multiwfn's default screening (RDG set to 100 where ρ ≥ 0.05 a.u.) would remove this very short hydrogen bond from the map, so the threshold was raised to 0.1 a.u.; covalent bond regions remain screened.\n\n")
hb_names = list(case.get('hbonds') or {})
hbond_desc = ('the blue disc' + ('s' if len(hb_names) > 1 else '') + ' at ' +
              ', '.join(h.replace('...', '···') for h in hb_names) +
              ': strong, attractive hydrogen bond' + ('s' if len(hb_names) > 1 else '')) if hb_names else 'the blue disc between the hydrogen-bond donor and acceptor'
L.append("*Draft interpretation.* The RDG scatter plot shows a spike at strongly negative sign(λ₂)ρ, which corresponds in the 3-D map to " + hbond_desc + ". "
         "Green, low-density regions between the two partners and around the rings mark van der Waals contacts, "
         "while the red discs at the centres of the aromatic rings indicate steric (ring-closure) repulsion. Only the hydrogen bonds inside this cluster are present; "
         "the remaining hydrogen bonds of the crystal involve further symmetry-related neighbours and cannot appear here." + hbtxt + "\n")
# ---- 5/6 ELF, LOL
for func, title, txt in [('elf', 'Electron localisation function (ELF)',
                          "*Draft interpretation.* ELF values close to 1 (red) appear on the X–H bonds and on the lone-pair regions of the electronegative atoms (oxygen, nitrogen, chlorine), showing highly localised electron pairs; "
                          "intermediate values (green) along the ring C–N and C–C bonds reflect delocalised π bonding within the rings, and values near zero (blue) at the ring centres and in the region between the fragments indicate electron depletion."),
                         ('lol', 'Localised orbital locator (LOL)',
                          "*Draft interpretation.* The LOL maps mirror the ELF picture with sharper contrast: high LOL (red) is concentrated on the covalent bond axes and the lone pairs, the delocalised π framework appears as continuous intermediate-LOL regions around the ring perimeters, and the low-LOL regions at the ring centres and in the hydrogen-bond region between the fragments show where electrons are delocalised or depleted.")]:
    L.append(f"## {5 if func == 'elf' else 6}. {title}\n")
    for plane in case.get('planes', {}):
        L.append(fig(f'{func}_{plane}.png', f'{func.upper()} colour-filled map in the {plane.replace("_", " ")} plane'))
    L.append(txt + '\n')
# ---- files
L.append("## Files\n\n`opt.xyz` optimised geometry · `geometry_vs_start.csv` bond lengths vs X-ray · `freq.json` frequencies · `fmo.json`/`fmo_table.md` descriptors · "
         "`mep_extrema.csv`/`mep_surface.json` ESP surface analysis · `nci_spikes.csv` RDG spikes · `fig/` all figures · `STAGES.md` stage log · `logs/` raw program output.\n")
(d / 'RESULTS.md').write_text(''.join(L)); stage_note(case, 'RESULTS.md assembled'); print('[report] wrote', d / 'RESULTS.md')
