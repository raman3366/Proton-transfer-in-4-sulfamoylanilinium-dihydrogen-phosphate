"""
2-D figures with matplotlib: MO energy-level diagram, RDG scatter (NCI), ELF and LOL colour-filled maps,
and a labelled structure diagram.  Usage: python dft/plots.py <case>
Inputs: results/<case>/{fmo.json, opt.xyz, rdg_scatter.txt, elf_<plane>.npz, lol_<plane>.npz}
"""
import sys, json, pathlib, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent)); from common import *

case = load_case(sys.argv[1]); D = case['dir']; FIG = D / 'fig'; FIG.mkdir(exist_ok=True)
plt.rcParams.update({'font.size': 11, 'font.family': 'DejaVu Sans', 'figure.dpi': 150, 'savefig.dpi': 300})

# 1. structure with atom labels (numbering used in all tables)
syms, xyz = read_xyz(D / 'opt.xyz'); lab = case.get('cif_labels') or atom_labels(syms)
c = xyz.mean(0); w, v = np.linalg.eigh(np.cov((xyz - c).T)); P = (xyz - c) @ v[:, [2, 1]]
col = {'C': '#404040', 'N': '#2050d0', 'O': '#d02020', 'H': '#b0b0b0'}
fig, ax = plt.subplots(figsize=(6, 5))
for i, j, d in bonds_from_geometry(syms, xyz): ax.plot(P[[i, j], 0], P[[i, j], 1], '-', color='#808080', lw=2, zorder=1)
for i, s in enumerate(syms):
    ax.scatter(P[i, 0], P[i, 1], s=350 if s != 'H' else 150, color=col.get(s, 'k'), zorder=2, edgecolor='k')
    ax.annotate(lab[i], P[i], xytext=(6, 6), textcoords='offset points', fontsize=8)
ax.set_aspect('equal'); ax.axis('off'); ax.set_title(f'Optimised geometry: {case["name"]} (B3LYP/6-311++G(d,p))', fontsize=10)
fig.tight_layout(); fig.savefig(FIG / 'structure_labels.png'); plt.close(fig)

# 2. MO energy-level diagram
if (D / 'fmo.json').exists():
  f = json.load(open(D / 'fmo.json')); e = np.array(f['mo_energies_eV']); n = f['n_occ']
  fig, ax = plt.subplots(figsize=(4.2, 5.5))
  for k in range(max(0, n - 6), min(len(e), n + 6)):
      occ = k < n; ax.hlines(e[k], 0.2, 0.8, color='#d02020' if occ else '#2050d0', lw=2.2)
      if k in (n - 1, n): ax.text(0.85, e[k], f'{"HOMO" if occ else "LUMO"}  {e[k]:.2f} eV', va='center', fontsize=9)
  ax.annotate('', xy=(0.5, e[n]), xytext=(0.5, e[n - 1]), arrowprops=dict(arrowstyle='<->', color='k'))
  ax.text(0.52, (e[n] + e[n - 1]) / 2, f'ΔE = {f["gap_eV"]:.2f} eV', fontsize=9, va='center')
  ax.set_xlim(0, 2.2); ax.set_xticks([]); ax.set_ylabel('Orbital energy (eV)'); ax.set_title('Frontier molecular orbitals', fontsize=10)
  fig.tight_layout(); fig.savefig(FIG / 'mo_levels.png'); plt.close(fig)

# 3. RDG scatter (NCI)
bgr = LinearSegmentedColormap.from_list('bgr', [(0, 0, 1), (0, 1, 0), (1, 0, 0)])
if (D / 'rdg_scatter.npz').exists():
    z = np.load(D / 'rdg_scatter.npz'); x0, y0 = z['x'], z['y']
    strong = (x0 < -0.05) & (y0 < 0.5)
    for name, xlim in [('rdg_scatter', (-0.05, 0.05))] + ([('rdg_scatter_wide', (float(np.floor(x0[strong].min() * 20) / 20) - 0.05, 0.05))] if strong.any() else []):
        m = (x0 > xlim[0]) & (x0 < xlim[1]) & (y0 < 2.0); x, y = x0[m], y0[m]
        if len(x) > 400000: idx = np.random.default_rng(0).choice(len(x), 400000, replace=False); x, y = x[idx], y[idx]
        fig, ax = plt.subplots(figsize=(6, 4.5))
        ax.scatter(x, y, c=x, cmap=bgr, vmin=-0.035, vmax=0.02, s=1.2, rasterized=True, linewidths=0)
        ax.set_xlim(*xlim); ax.set_ylim(0, 2); ax.set_xlabel(r'sign($\lambda_2$)$\rho$ (a.u.)'); ax.set_ylabel('RDG (a.u.)')
        ax.axhline(0.5, color='k', ls='--', lw=0.8); ax.set_title('Reduced density gradient vs sign(λ2)ρ' + (' (extended range)' if 'wide' in name else ''), fontsize=10)
        sm = plt.cm.ScalarMappable(cmap=bgr, norm=plt.Normalize(-0.035, 0.02)); fig.colorbar(sm, ax=ax, label='sign(λ2)ρ  (blue: H-bond, green: vdW, red: steric)')
        fig.tight_layout(); fig.savefig(FIG / f'{name}.png'); plt.close(fig)

# 4. ELF / LOL colour-filled contour maps for each requested plane
for func in ('elf', 'lol'):
    for plane in case.get('planes', {}):
        p = D / f'{func}_{plane}.npz'
        if not p.exists(): continue
        z = np.load(p); X, Y, F = z['X'], z['Y'], z['F']
        fig, ax = plt.subplots(figsize=(6.2, 5))
        cf = ax.contourf(X, Y, F, levels=np.linspace(0, 1, 41), cmap='jet'); ax.contour(X, Y, F, levels=np.linspace(0.1, 0.9, 9), colors='k', linewidths=0.3)
        for i, s in enumerate(syms):
            if 'atoms_xy' in z and z['in_plane'][i]: ax.text(z['atoms_xy'][i, 0], z['atoms_xy'][i, 1], lab[i], fontsize=8, ha='center', va='center', color='w', fontweight='bold')
        fig.colorbar(cf, ax=ax, ticks=np.linspace(0, 1, 6), label=func.upper()); ax.set_aspect('equal')
        ax.set_xlabel('x (Å)'); ax.set_ylabel('y (Å)'); ax.set_title(f'{func.upper()} in the {plane.replace("_", " ")} plane', fontsize=10)
        fig.tight_layout(); fig.savefig(FIG / f'{func}_{plane}.png'); plt.close(fig)
print('[plots] figures written to', FIG)
