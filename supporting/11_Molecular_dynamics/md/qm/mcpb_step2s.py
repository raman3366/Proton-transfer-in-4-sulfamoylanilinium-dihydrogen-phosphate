"""MCPB.py step 2s (Seminario) with a guard against a fault of MCPB.py 2026 under NumPy 2.

    AMBERHOME=/root/cqc-sfa-tools/md /root/cqc-sfa-tools/md/bin/python md/qm/mcpb_step2s.py CA2.in   (in the work folder)

The fault.  Seminario's method diagonalises the 3 x 3 block of the Hessian between two atoms.  The block is
not symmetric, so two of its eigenvalues can be a complex-conjugate pair.  MCPB.py sums
eigenvalue * |direction . eigenvector| over the three; for a conjugate pair the imaginary parts cancel, but
the sum keeps a complex type, and `round()` of a NumPy 2 complex number raises
    TypeError: type numpy.complex128 doesn't define __round__ method
(pymsmt/mcpb/gene_final_frcmod_file.py, avg_bond_para and avg_angle_para).  Found with the synthetic test
md/qm/test_after_qm.py; it will also happen with real Hessians whenever a block has such a pair.

The guard.  The two averaging functions are wrapped: every force constant and equilibrium value is replaced
by its real part before MCPB.py rounds it, after checking that the imaginary part is negligible (below
1e-6 of the modulus; otherwise this script stops).  Nothing else of MCPB.py is changed, and where no complex
number occurs the result is identical to plain `MCPB.py -s 2s`.  What was converted is printed.
"""
import runpy
import sys

import numpy

import pymsmt.mcpb.gene_final_frcmod_file as g

_bond, _angle = g.avg_bond_para, g.avg_angle_para
_seen = dict(n=0, complex=0, largest_relative_imaginary=0.0)


def _real(values, label):
    out = []
    for v in values:
        _seen['n'] += 1
        if numpy.iscomplexobj(v):
            _seen['complex'] += 1
            rel = abs(numpy.imag(v)) / max(abs(v), 1e-300)
            _seen['largest_relative_imaginary'] = max(_seen['largest_relative_imaginary'], float(rel))
            if rel > 1e-6:
                raise SystemExit(f'mcpb_step2s.py: {label} has a real imaginary part ({v}); not a rounding effect')
            v = float(numpy.real(v))
        out.append(v)
    return out


def avg_bond_para(misbondat12, bondlen, bfconst, sigf=4):
    return _bond(misbondat12, _real(bondlen, misbondat12), _real(bfconst, misbondat12), sigf)


def avg_angle_para(misangat123, angvals, afconst):
    return _angle(misangat123, _real(angvals, misangat123), _real(afconst, misangat123))


g.avg_bond_para, g.avg_angle_para = avg_bond_para, avg_angle_para

if __name__ == '__main__':
    inp = sys.argv[1]
    sys.argv = ['MCPB.py', '-i', inp, '-s', '2s']
    try:
        runpy.run_path(sys.prefix + '/bin/MCPB.py', run_name='__main__')
    except SystemExit as e:
        if e.code not in (None, 0):
            raise
    print(f"mcpb_step2s.py: {_seen['n']} values passed to MCPB.py's averaging; {_seen['complex']} were of complex "
          f"type, largest |imaginary|/|value| {_seen['largest_relative_imaginary']:.1e} (real part taken)")
