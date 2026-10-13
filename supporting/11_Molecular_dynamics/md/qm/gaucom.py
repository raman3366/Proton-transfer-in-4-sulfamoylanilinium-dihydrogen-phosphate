"""Adapter (a), first half: read a Gaussian input file (.com) as MCPB.py writes it.

Pure Python, no third-party import, so that it runs in the quantum-chemistry environment and in the MD
environment alike.  Only the constructs MCPB.py 2026 writes are understood; anything else stops with an
error instead of being guessed:

    %Chk=..., %Mem=..., %NProcShared=...          link-0 lines (kept, not used)
    # B3LYP/6-31G* ... [second route line]        route; the method and the basis are read from it
    blank, title (CLR), blank
    charge  multiplicity
    El  x y z            (small model, large model with large_opt 0 or 2)
    El  flag x y z       (large model with large_opt 1: flag 0 = free, -1 = frozen)
    blank
    El radius !comment   (after Pop(MK,ReadRadii): radius of the metal ion for the fit points, in Angstrom)

The force-constant input (<group>_small_fc.com) has no geometry (Geom=AllCheckpoint): `atoms` is then empty
and `geometry_from_checkpoint` is True; the caller takes the optimised geometry of the preceding job.
"""
import pathlib
import re

BOHR = 0.52917721092  # Angstrom per Bohr, CODATA 2010 (the value PySCF 2.14 uses: pyscf.lib.param.BOHR)


class ComFile:
    def __init__(self):
        self.path = None
        self.link0 = {}
        self.route = ''
        self.method = None
        self.basis = None
        self.title = ''
        self.charge = None
        self.multiplicity = None
        self.symbols = []
        self.coords = []          # Angstrom
        self.frozen = []          # True where the optimisation flag is -1
        self.radii = {}           # element symbol -> radius in Angstrom (ReadRadii)
        self.geometry_from_checkpoint = False
        self.keywords = {}

    @property
    def natm(self):
        return len(self.symbols)


def _keywords(route):
    kw = {}
    r = route.lower()
    kw['opt'] = bool(re.search(r'(^|\s)opt(\s|=|\(|$)', r))
    kw['opt_calcall'] = 'opt=calcall' in r or 'calcall' in r
    kw['pop_mk'] = bool(re.search(r'pop\(?=?\(?mk', r))
    kw['read_radii'] = 'readradii' in r
    kw['grid'] = 'ultrafine' if 'ultrafine' in r else None
    m = re.search(r'iop\(([^)]*)\)', r)
    kw['iop'] = {}
    for m in re.finditer(r'iop\(([^)]*)\)', r):
        for item in m.group(1).split(','):
            k, v = item.split('=')
            kw['iop'][k.strip()] = int(v)
    kw['geom_allcheckpoint'] = 'geom=allcheckpoint' in r
    return kw


def read_com(path):
    path = pathlib.Path(path)
    lines = path.read_text().splitlines()
    c = ComFile()
    c.path = str(path)
    i = 0
    while i < len(lines) and lines[i].startswith('%'):
        k, _, v = lines[i][1:].partition('=')
        c.link0[k.strip().lower()] = v.strip()
        i += 1
    if not lines[i].lstrip().startswith('#'):
        raise ValueError(f'{path}: route line expected at line {i + 1}')
    route = []
    while i < len(lines) and lines[i].strip():
        route.append(lines[i].strip())
        i += 1
    c.route = ' '.join(route)
    m = re.match(r'#\s*([A-Za-z0-9]+)/(\S+)', c.route)
    if not m:
        raise ValueError(f'{path}: cannot read method/basis from the route "{c.route}"')
    c.method, c.basis = m.group(1), m.group(2)
    c.keywords = _keywords(c.route)
    if c.keywords['geom_allcheckpoint']:
        c.geometry_from_checkpoint = True
        return c
    i += 1                      # blank
    c.title = lines[i].strip()
    i += 2                      # title, blank
    t = lines[i].split()
    if len(t) != 2:
        raise ValueError(f'{path}: charge and multiplicity expected at line {i + 1}')
    c.charge, c.multiplicity = int(t[0]), int(t[1])
    i += 1
    while i < len(lines) and lines[i].strip():
        t = lines[i].split()
        if len(t) == 4:
            flag = 0
        elif len(t) == 5:
            flag = int(t[1])
            if flag not in (0, -1):
                raise ValueError(f'{path}: optimisation flag {flag} at line {i + 1}')
        else:
            raise ValueError(f'{path}: atom line not understood at line {i + 1}: {lines[i]!r}')
        sym = t[0][0].upper() + t[0][1:].lower()
        c.symbols.append(sym)
        c.coords.append(tuple(float(x) for x in t[-3:]))
        c.frozen.append(flag == -1)
        i += 1
    # what follows the geometry: radii (ReadRadii) or ModRedundant lines (not written with our settings)
    for l in lines[i:]:
        if not l.strip():
            continue
        t = l.split('!')[0].split()
        if len(t) == 2 and re.fullmatch(r'[A-Za-z]{1,2}', t[0]):
            c.radii[t[0][0].upper() + t[0][1:].lower()] = float(t[1])
        else:
            raise ValueError(f'{path}: line after the geometry not understood: {l!r}')
    if c.keywords['read_radii'] and not c.radii:
        raise ValueError(f'{path}: Pop(MK,ReadRadii) but no radius found')
    return c


if __name__ == '__main__':
    import json
    import sys
    for f in sys.argv[1:]:
        c = read_com(f)
        comp = {}
        for s in c.symbols:
            comp[s] = comp.get(s, 0) + 1
        print(json.dumps(dict(file=f, method=c.method, basis=c.basis, charge=c.charge, mult=c.multiplicity,
                              natm=c.natm, composition=comp, n_frozen=sum(c.frozen), radii=c.radii,
                              keywords=c.keywords, from_chk=c.geometry_from_checkpoint)))
