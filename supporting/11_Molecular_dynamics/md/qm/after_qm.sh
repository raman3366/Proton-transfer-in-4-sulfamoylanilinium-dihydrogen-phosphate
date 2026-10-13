#!/bin/bash
# Stage S2, after the quantum-chemical jobs: MCPB.py steps 2, 3 and 4, the tleap run and the checks.
#
#   md/qm/after_qm.sh [step2] [workdir] [reportdir]
#
#   step2     2s  Seminario from <group>_small_opt.fchk (the route of the plan; default)
#             2e  empirical zinc constants (Yu, Li and Merz 2018), fallback (a) of md/PLAN.md, needs no fchk
#   workdir   the MCPB.py folder (default: md/qm/mcpb)
#   reportdir where the check writes its JSON (default: md/qm/results)
#
# Needs in the work folder: the files of step 1 (md/qm/mcpb_step1.py), CA2_large_mk.log (always) and
# CA2_small_opt.fchk (for 2s), both written by md/qm/run_qm.py.  Nothing here is random.  Every step writes
# its screen output to a log in the work folder and the script stops at the first error.
#
# HOW THE REPAIRED LIGAND TERMS ARRIVE.  CA2.in names `gaff 2` and
#   frcmod_files sfa_anion_gaff2.frcmod sfa_anion_repair.frcmod      (in this order)
# MCPB.py step 2 renames the type of the ligating sulfonamide nitrogen (n2 -> Y4) and copies every bond,
# angle and torsion that carries it from gaff2.dat and these two files, the later file replacing the
# earlier, into CA2_mcpbpy.frcmod under the new type (Y4-sy, ca-sy-Y4, hn-Y4-sy, Y4-sy-o, ca-sy-Y4-hn,
# hn-Y4-sy-o, Y4-hn).  tleap then loads the two ligand files (terms without Y4: the ring, the aniline) and
# CA2_mcpbpy.frcmod (terms with Y4, and the zinc).  md/qm/after_qm_check.py compares what arrived with the
# table in md/LIGAND.md, section 6, and stops this script if one value differs.
#
# CHARGES.  After MCPB.py step 3 the fit is repeated with the six pairs of chemically equal ligand atoms made
# equal (md/qm/resp_equiv.py; added on 6 Oct 2026 for gate G-R, md/GR.md).  The mol2 files that the later
# steps and the system builder read carry these constrained charges; MCPB.py's own are in
# charges_unconstrained/.
#
# THREE FAULTS OF MCPB.py 2026 THAT THIS SCRIPT WORKS AROUND
#  1. With `gaff 2` step 4 prints "source leaprc.gaff2" to the screen instead of into the tleap input
#     (pymsmt/mcpb/amber_modeling.py, line 240).  The line is added here.
#  2. Step 4 sources oldff/leaprc.ff14SB.  The system of this project is built with
#     leaprc.protein.ff14SB (md/check_tleap.py); the input written here uses that file, so that the
#     protein has the same templates in both.
#  3. Step 2s stops with "TypeError: type numpy.complex128 doesn't define __round__ method" when a 3 x 3
#     block of the Hessian has a complex-conjugate pair of eigenvalues (NumPy 2).  Step 2s is therefore run
#     through md/qm/mcpb_step2s.py, which takes the real part (the imaginary parts cancel in the sum) and
#     prints how many values it touched.
set -euo pipefail

STEP2=${1:-2s}
HERE=$(cd "$(dirname "$0")" && pwd)
WORK=${2:-$HERE/mcpb}
REPORT=${3:-$HERE/results}
GROUP=CA2
export AMBERHOME=/root/cqc-sfa-tools/md
export PATH=$AMBERHOME/bin:$PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=$AMBERHOME/bin/python

case "$STEP2" in 2s|2e) ;; *) echo "step2 must be 2s or 2e" >&2; exit 2;; esac
cd "$WORK"
mkdir -p "$REPORT"
test -s ${GROUP}_large_mk.log || { echo "missing ${GROUP}_large_mk.log (run md/qm/run_qm.py, job esp)" >&2; exit 1; }
if [ "$STEP2" = 2s ]; then
  test -s ${GROUP}_small_opt.fchk || { echo "missing ${GROUP}_small_opt.fchk (jobs opt and hess); or use 2e" >&2; exit 1; }
fi

echo "== MCPB.py step $STEP2: bonded terms"
rm -f ${GROUP}_mcpbpy.frcmod ${GROUP}_mcpbpy_pre.frcmod
if [ "$STEP2" = 2s ]; then
  # through md/qm/mcpb_step2s.py: MCPB.py 2026 stops with a TypeError under NumPy 2 when a block of the
  # Hessian has a complex pair of eigenvalues (fault 3 below the header); the wrapper takes the real part
  nice -n 10 $PY "$HERE/mcpb_step2s.py" ${GROUP}.in > step2s.log 2>&1 || { tail -5 step2s.log >&2; exit 1; }
else
  nice -n 10 $PY $AMBERHOME/bin/MCPB.py -i ${GROUP}.in -s 2e > step2e.log 2>&1 || { tail -5 step2e.log >&2; exit 1; }
fi
test -s ${GROUP}_mcpbpy.frcmod

echo "== MCPB.py step 3 (3b): RESP charges of the five residues of the site"
# espgen step inside MCPB.py reads ${GROUP}_large_mk.log and writes ${GROUP}_large_mk.esp; resp runs twice
rm -f HD1.mol2 HD2.mol2 HE1.mol2 ZN1.mol2 SN1.mol2
nice -n 10 $PY $AMBERHOME/bin/MCPB.py -i ${GROUP}.in -s 3 > step3.log 2>&1 || { tail -5 step3.log >&2; exit 1; }
for r in HD1 HD2 HE1 ZN1 SN1; do test -s $r.mol2 || { echo "step 3 did not write $r.mol2" >&2; exit 1; }; done

echo "== RESP again with the chemically equal ligand atoms made equal (md/PLAN.md, 6 Oct 2026, 3:02 AM, point 5)"
# md/qm/resp_equiv.py: O1=O2, C2=C6, C3=C5, H2=H6, H3=H5, HN11=HN12 in the first resp stage.  It rewrites the
# five mol2 files with the constrained charges (USED from here on) and keeps the files MCPB.py wrote in
# charges_unconstrained/ (reported beside).  Report: resp_equiv.json here and in the report folder.
nice -n 10 $PY "$HERE/resp_equiv.py" "$WORK" "$REPORT/resp_equiv_${STEP2}.json" > resp_equiv.log 2>&1 || { tail -8 resp_equiv.log >&2; exit 1; }
tail -6 resp_equiv.log

echo "== MCPB.py step 4: tleap input and the PDB file with the renamed residues"
rm -f ${GROUP}_tleap.in ${GROUP}_mcpbpy.pdb
nice -n 10 $PY $AMBERHOME/bin/MCPB.py -i ${GROUP}.in -s 4 > step4.log 2>&1 || { tail -5 step4.log >&2; exit 1; }
test -s ${GROUP}_tleap.in && test -s ${GROUP}_mcpbpy.pdb

echo "== tleap input for the dry complex of start X (protein, zinc, ligand; no water, no ions)"
# from the input of MCPB.py: other protein leaprc, GAFF2 added, everything after the first saveamberparm cut
$PY - <<EOF
import re, pathlib
src = pathlib.Path('${GROUP}_tleap.in').read_text().splitlines()
out, done = [], False
for l in src:
    if l.startswith('source oldff/leaprc.ff14SB'):
        out += ['source leaprc.protein.ff14SB', 'source leaprc.gaff2']
        continue
    if l.strip() == 'source leaprc.gaff2':
        continue                                   # should a later MCPB.py write it itself
    if l.startswith(('solvatebox', 'addions', 'quit')) or '_solv' in l:
        continue
    out.append(l)
assert any(l.startswith('source leaprc.gaff2') for l in out)
assert sum(l.startswith('bond ') for l in out) == 10, 'expected 4 bonds to zinc and 6 peptide bonds'
out += ['check mol', 'charge mol', 'quit']
pathlib.Path('${GROUP}_tleap_dry.in').write_text('\n'.join(out) + '\n')
EOF
rm -f ${GROUP}_dry.prmtop ${GROUP}_dry.inpcrd ${GROUP}_dry.pdb
nice -n 10 tleap -s -f ${GROUP}_tleap_dry.in > ${GROUP}_tleap_dry.log 2>&1 || true
test -s ${GROUP}_dry.prmtop || { echo "tleap wrote no topology, see $WORK/${GROUP}_tleap_dry.log" >&2; exit 1; }

echo "== checks (terms carried over, zinc terms, charges, MCPB.py's own criteria, a short minimisation)"
nice -n 10 $PY "$HERE/after_qm_check.py" "$WORK" "$REPORT" "$STEP2"
echo "== done: $REPORT/after_qm_check_${STEP2}.json"
echo "   for the system builder (stage S3): HD1.mol2 HD2.mol2 HE1.mol2 ZN1.mol2 SN1.mol2, ${GROUP}_mcpbpy.frcmod,"
echo "   the two ligand frcmod files, the addAtomTypes block and the ten bond commands of ${GROUP}_tleap_dry.in."
