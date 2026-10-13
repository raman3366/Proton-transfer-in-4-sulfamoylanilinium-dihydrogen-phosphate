"""
3-D isosurface renders with ParaView (pvbatch), replacing VMD.
Cube files are converted here to legacy VTK image data in Angstrom (ParaView's own cube reader drops origin/spacing),
and the molecule is drawn explicitly as spheres + tubes from opt.xyz.
Usage: pvbatch --force-offscreen-rendering dft/render.py <case>     (spawns one pvbatch per figure)
Produces results/<case>/fig/: homo.png, lumo.png, mep.png, nci.png (the latter two if the input cubes exist)
"""
import sys, os, json, pathlib, subprocess
import numpy as np

case = sys.argv[-1]   # case name is always the last argument
ROOT = pathlib.Path(__file__).resolve().parent.parent
D = next((p.parent.parent / 'results' / case for p in (ROOT/'molecules').glob(f'*/cases/{case}.json')),
         ROOT / 'results' / case); FIG = D / 'fig'; FIG.mkdir(exist_ok=True); VTK = D / 'vtk'; VTK.mkdir(exist_ok=True)
RES = [1600, 1200]; BOHR = 0.529177210903
PVBATCH = os.environ.get('PVBATCH') or (str(ROOT / 'tools' / 'bin' / 'pvbatch')
           if (ROOT / 'tools' / 'bin' / 'pvbatch').exists()
           else '/Applications/ParaView-6.1.1.app/Contents/bin/pvbatch')
COLORS = {'C': (0.30, 0.30, 0.30), 'N': (0.15, 0.30, 0.85), 'O': (0.85, 0.15, 0.15), 'H': (0.85, 0.85, 0.85),
          'S': (0.90, 0.80, 0.20), 'Cl': (0.12, 0.80, 0.12), 'Br': (0.60, 0.20, 0.10), 'F': (0.50, 0.85, 0.85),
          'P': (1.00, 0.55, 0.10), 'Zn': (0.55, 0.50, 0.65), 'Cu': (0.78, 0.50, 0.20), 'Fe': (0.88, 0.40, 0.20),
          'Ni': (0.31, 0.82, 0.31), 'Co': (0.94, 0.56, 0.63), 'Mn': (0.61, 0.48, 0.78), 'Cd': (1.00, 0.85, 0.56),
          'Na': (0.67, 0.36, 0.95), 'K': (0.56, 0.25, 0.83), 'Mg': (0.54, 1.00, 0.00), 'Ca': (0.24, 1.00, 0.00)}
DEFAULT_COLOR = (0.60, 0.60, 0.60)
RADII = {'C': 0.36, 'N': 0.36, 'O': 0.36, 'H': 0.22, 'S': 0.42, 'Cl': 0.40, 'Br': 0.44, 'F': 0.32,
         'P': 0.42, 'Zn': 0.46, 'Cu': 0.46, 'Fe': 0.48, 'Ni': 0.46, 'Co': 0.46, 'Mn': 0.48, 'Cd': 0.50,
         'Na': 0.48, 'K': 0.56, 'Mg': 0.44, 'Ca': 0.52}
DEFAULT_RADIUS = 0.38
COV = {'H': 0.31, 'B': 0.84, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'Na': 1.66, 'Mg': 1.41,
       'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'K': 2.03, 'Ca': 1.76, 'Mn': 1.39,
       'Fe': 1.32, 'Co': 1.26, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 'Br': 1.20, 'Cd': 1.44, 'I': 1.39}

def read_xyz(p):
    L = open(p).read().strip().splitlines(); n = int(L[0])
    return [l.split()[0] for l in L[2:2+n]], np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2+n]])

def cube_to_vtk(cube, out):
    """Gaussian cube (Bohr) -> legacy VTK STRUCTURED_POINTS in Angstrom, x fastest as VTK requires."""
    if out.exists() and out.stat().st_mtime > cube.stat().st_mtime: return out
    L = open(cube).read().split('\n'); nat = int(L[2].split()[0]); org = np.array([float(v) for v in L[2].split()[1:4]])
    n = []; ax = []
    for i in (3, 4, 5):
        p = L[i].split(); n.append(int(p[0])); ax.append([float(v) for v in p[1:4]])
    ax = np.array(ax); skip = 6 + abs(nat) + (1 if nat < 0 else 0)
    vals = np.array(' '.join(L[skip:]).split(), float)[: n[0]*n[1]*n[2]].reshape(n[0], n[1], n[2])   # cube: z fastest
    sp = np.array([ax[0, 0], ax[1, 1], ax[2, 2]]); org_A, sp_A = org * BOHR, sp * BOHR
    if np.any(sp <= 0): raise ValueError('non-positive cube spacing')
    with open(out, 'w') as f:
        f.write(f'# vtk DataFile Version 3.0\n{cube.name}\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS {n[0]} {n[1]} {n[2]}\n'
                f'ORIGIN {org_A[0]:.6f} {org_A[1]:.6f} {org_A[2]:.6f}\nSPACING {sp_A[0]:.6f} {sp_A[1]:.6f} {sp_A[2]:.6f}\n'
                f'POINT_DATA {vals.size}\nSCALARS value float 1\nLOOKUP_TABLE default\n')
        np.savetxt(f, vals.transpose(2, 1, 0).ravel(), fmt='%.6e')       # x fastest
    return out

def molecule_vtk():
    """Per-element atom point sets (with radius) and a bond line set, as legacy VTK polydata."""
    syms, xyz = read_xyz(D / 'opt.xyz'); files = {}
    for el in sorted(set(syms)):
        idx = [i for i, s in enumerate(syms) if s == el]; p = VTK / f'atoms_{el}.vtk'
        with open(p, 'w') as f:
            f.write(f'# vtk DataFile Version 3.0\natoms {el}\nASCII\nDATASET POLYDATA\nPOINTS {len(idx)} float\n')
            for i in idx: f.write(f'{xyz[i,0]:.6f} {xyz[i,1]:.6f} {xyz[i,2]:.6f}\n')
            f.write(f'VERTICES {len(idx)} {2*len(idx)}\n' + ''.join(f'1 {k}\n' for k in range(len(idx))))
        files[el] = p
    bonds = [(i, j) for i in range(len(syms)) for j in range(i+1, len(syms))
             if np.linalg.norm(xyz[i]-xyz[j]) < 1.25 * (COV.get(syms[i], .8) + COV.get(syms[j], .8))]
    p = VTK / 'bonds.vtk'
    with open(p, 'w') as f:
        f.write(f'# vtk DataFile Version 3.0\nbonds\nASCII\nDATASET POLYDATA\nPOINTS {len(syms)} float\n')
        for x in xyz: f.write(f'{x[0]:.6f} {x[1]:.6f} {x[2]:.6f}\n')
        f.write(f'LINES {len(bonds)} {3*len(bonds)}\n' + ''.join(f'2 {i} {j}\n' for i, j in bonds))
    return files, p, xyz

# ---------------- ParaView part (runs inside pvbatch) ----------------
def pv_setup():
    from paraview.simple import GetActiveViewOrCreate
    v = GetActiveViewOrCreate('RenderView'); v.ViewSize = RES; v.OrientationAxesVisibility = 0; v.Background = [1, 1, 1]
    try: v.UseColorPaletteForBackground = 0
    except Exception: pass
    return v

def pv_molecule(view):
    from paraview.simple import LegacyVTKReader, Glyph, Tube, Show, ColorBy
    files, bonds, xyz = molecule_vtk()
    for el, p in files.items():
        r = LegacyVTKReader(FileNames=[str(p)])
        g = Glyph(Input=r, GlyphType='Sphere', GlyphMode='All Points', ScaleArray=['POINTS', 'No scale array'], ScaleFactor=2 * RADII.get(el, 0.3))
        g.GlyphType.ThetaResolution = 32; g.GlyphType.PhiResolution = 32
        rep = Show(g, view); rep.ColorArrayName = ['POINTS', '']; rep.DiffuseColor = list(COLORS.get(el, (0.5, 0.5, 0.5))); rep.AmbientColor = rep.DiffuseColor
        rep.Specular = 0.3
    t = Tube(Input=LegacyVTKReader(FileNames=[str(bonds)]), Radius=0.10, NumberofSides=24)
    rep = Show(t, view); rep.ColorArrayName = ['POINTS', '']; rep.DiffuseColor = [0.55, 0.55, 0.55]; rep.AmbientColor = rep.DiffuseColor
    return xyz

def pv_finish(view, xyz, out, zoom=1.0):
    from paraview.simple import Render, SaveScreenshot
    c = xyz.mean(0); w, v = np.linalg.eigh(np.cov((xyz - c).T)); normal, up = v[:, 0], v[:, 2]
    if normal[2] < 0: normal = -normal
    ext = np.ptp(xyz, axis=0).max() + 5
    view.CameraFocalPoint = list(c); view.CameraPosition = list(c + normal * ext * 3); view.CameraViewUp = list(up)
    view.CameraParallelProjection = 1; view.CameraParallelScale = ext * 0.42 / zoom
    Render(view); SaveScreenshot(str(out), view, ImageResolution=RES, TransparentBackground=0); print('wrote', out, flush=True)

def orbital_png(cube, out, iso=0.02):
    from paraview.simple import LegacyVTKReader, Contour, Show, ColorBy
    view = pv_setup(); xyz = pv_molecule(view)
    r = LegacyVTKReader(FileNames=[str(cube_to_vtk(cube, VTK / (cube.stem + '.vtk')))])
    for val, col in [(iso, [0.10, 0.35, 0.95]), (-iso, [0.90, 0.15, 0.15])]:
        ct = Contour(Input=r, ContourBy=['POINTS', 'value'], Isosurfaces=[val])
        rep = Show(ct, view); rep.ColorArrayName = ['POINTS', '']; rep.DiffuseColor = col; rep.AmbientColor = col; rep.Opacity = 0.6
    pv_finish(view, xyz, out)

def mapped_png(iso_cube, color_cube, iso, out, rng, rgb_points, title, opacity=1.0, molecule=True):
    from paraview.simple import LegacyVTKReader, Contour, Calculator, ResampleWithDataset, Show, ColorBy, GetColorTransferFunction, GetScalarBar
    view = pv_setup(); xyz = pv_molecule(view) if molecule else read_xyz(D / 'opt.xyz')[1]
    r1 = LegacyVTKReader(FileNames=[str(cube_to_vtk(iso_cube, VTK / (iso_cube.stem + '.vtk')))])
    r2 = LegacyVTKReader(FileNames=[str(cube_to_vtk(color_cube, VTK / (color_cube.stem + '.vtk')))])
    ct = Contour(Input=r1, ContourBy=['POINTS', 'value'], Isosurfaces=[iso])
    c2 = Calculator(Input=r2, ResultArrayName='colour', Function='value')
    coloured = ResampleWithDataset(SourceDataArrays=c2, DestinationMesh=ct); coloured.PassPointArrays = 1
    rep = Show(coloured, view); ColorBy(rep, ('POINTS', 'colour')); rep.Opacity = opacity; rep.Specular = 0.2
    lut = GetColorTransferFunction('colour'); lut.RGBPoints = rgb_points; lut.ColorSpace = 'RGB'; lut.RescaleTransferFunction(rng[0], rng[1])
    bar = GetScalarBar(lut, view); bar.Title = title; bar.ComponentTitle = ''; bar.TitleColor = [0, 0, 0]; bar.LabelColor = [0, 0, 0]
    bar.Visibility = 1; bar.WindowLocation = 'Lower Right Corner'; bar.ScalarBarLength = 0.33; bar.TitleFontSize = 20; bar.LabelFontSize = 18
    pv_finish(view, xyz, out)

def points(vals, cols): return [x for v, c in zip(vals, cols) for x in (v, *c)]
MEP_COLS = [(1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)]      # red (negative) ... blue (positive)
NCI_COLS = [(0, 0, 1), (0, 1, 0), (1, 0, 0)]                                # blue (H-bond) green (vdW) red (steric)

def one(kind):
    if kind == 'homo': orbital_png(D / 'homo.cube', FIG / 'homo.png')
    elif kind == 'lumo': orbital_png(D / 'lumo.cube', FIG / 'lumo.png')
    elif kind == 'mep':
        mep = json.load(open(D / 'mep_range.json')) if (D / 'mep_range.json').exists() else {'vmin': -0.05, 'vmax': 0.05}
        lo, hi = mep['vmin'], mep['vmax']      # actual surface range: red at the most negative, blue at the most positive, green at 0
        mapped_png(D / 'density.cube', D / 'esp.cube', 0.001, FIG / 'mep.png', (lo, hi), points([lo, lo / 2, 0.0, hi / 2, hi], MEP_COLS), 'ESP (a.u.)', opacity=0.75, molecule=True)
    elif kind == 'nci':
        mapped_png(D / 'rdg.cube', D / 'sl2r.cube', 0.5, FIG / 'nci.png', (-0.035, 0.02), points([-0.035, 0.0, 0.02], NCI_COLS), 'sign(l2)rho (a.u.)')

if __name__ == '__main__':
    if '--single' in sys.argv:
        one(sys.argv[sys.argv.index('--single') + 1])
    else:
        wanted = [('homo', ['homo.cube']), ('lumo', ['lumo.cube']), ('mep', ['density.cube', 'esp.cube']), ('nci', ['rdg.cube', 'sl2r.cube'])]
        for kind, needs in wanted:
            if all((D / f).exists() for f in needs):
                rc = subprocess.run([PVBATCH, '--force-offscreen-rendering', __file__, '--single', kind, case]).returncode
                print(f'[render] {kind}: {"ok" if rc == 0 else "FAILED"}', flush=True)
            else:
                print(f'[render] {kind}: skipped (missing {[f for f in needs if not (D / f).exists()]})')
