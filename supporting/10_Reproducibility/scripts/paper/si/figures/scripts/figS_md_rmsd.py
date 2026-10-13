"""SI figure (figS_md_rmsd): the ligand RMSD of the eleven molecular-dynamics runs.

    nice -n 10 tools/paper-env/bin/python paper/si/figures/scripts/figS_md_rmsd.py

(a-k) RMSD of the 11 heavy atoms of the ligand from the crystal pose of 6RL9 after the fit of the 258 C-alpha
      atoms (smallest of the four equivalent labellings), in the layout of figS_md_orientation: the nine runs with
      zinc set 2e in a 3 x 3 block (rows: start X, N, R) and the two sensitivity runs of start R below.  Grey
      strip: the equilibration after the release of the ligand, stretched.  Broken line: the 2.0 A of rule 1,
      which reads the mean of the last half of a crystal-start run.  Triangles as in figS_md_orientation.
Data: md/results/series/ (rule 6 of md/PLAN.md); verdicts and times: facts md.* (section 40).
Drawing code: figS_md_common.py.
"""
import figS_md_common as C

if __name__ == '__main__':
    C.make('rmsd', 'figS_md_rmsd')
