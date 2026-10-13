#!/usr/bin/env bash
# Gate G2: the full redocking / cross-docking matrix on the Pi.
#
# Redocking is the only evidence that a docking protocol reproduces a known answer, so it runs before
# any production docking of the DFT geometries.  Each job is 10 seeds at exhaustiveness 32 with its own
# starting conformer; jobs run one at a time with --cpu 3 and nice so the Pi stays usable.
#
# The four extra receptors are sensitivity cases, not duplicates:
#   T1_3TZF_seed12345  the 66-67 loop rebuilt from a different PDBFixer seed
#   T2_1AJ0_noSO4      the sulfate in the pyrophosphate subsite removed
set -u
cd "$(dirname "$0")/.."
PY=/root/cqc-sfa-tools/dock/bin/python
export PATH=/root/cqc-sfa-tools/dock/bin:$PATH

run() {  # target ligand scoring [extra...]
  echo "=== $1 $2 $3 ${4:-}"
  nice -n 10 "$PY" dock/run_vina.py --target "$1" --ligand "$2" --scoring "$3" \
      --seeds 10 --exhaustiveness 32 --cpu 3 ${4:-} 2>&1 | grep -v 'tagged as 2D'
}

# --- redocking: the co-crystal ligand back into its own structure -----------------------------
run T1_3TZF   smx_anion   vina
run T1_3TZF   smx_neutral vina
run T2_1AJ0   sfa_neutral vina
run T3_6RL9   sfa_neutral vina
run T3_6RL9   sfa_anion   vina
run T4_3HS4   azm_neutral vina
run T4_3HS4   azm_anion   vina
run T5_7S2K   dhp_anion   vina

# --- cross-docking: the other structure's ligand ----------------------------------------------
run T3_6RL9   azm_neutral vina
run T3_6RL9   azm_anion   vina
run T4_3HS4   sfa_neutral vina
run T4_3HS4   sfa_anion   vina

# --- receptor-preparation sensitivity ---------------------------------------------------------
run T1_3TZF_seed12345 smx_anion   vina
run T2_1AJ0_noSO4     sfa_neutral vina

# --- second scoring function on the primary redocks -------------------------------------------
run T1_3TZF   smx_anion   vinardo
run T2_1AJ0   sfa_neutral vinardo
run T3_6RL9   sfa_neutral vinardo
run T4_3HS4   azm_neutral vinardo
run T5_7S2K   dhp_anion   vinardo

# --- control: how much easier is it if the search starts from the crystal conformer? -----------
run T3_6RL9   sfa_neutral vina --from-crystal
echo "ALL DONE"
