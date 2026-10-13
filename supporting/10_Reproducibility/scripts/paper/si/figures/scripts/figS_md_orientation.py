"""SI figure (figS_md_orientation): the orientation coordinate of the eleven molecular-dynamics runs.

    nice -n 10 tools/paper-env/bin/python paper/si/figures/scripts/figS_md_orientation.py

(a-i) d1 (aniline N1 to Phe131 CZ) and d2 (N1 to Asn67 OD1) against time for the nine runs of 100 ns with zinc
      set 2e: rows start X (crystal pose), N (docked near-native pose) and R (docked ring-rotated pose), columns
      runs 1 to 3.  Grey strip: the 1.3 ns of the equilibration after the release of the ligand, stretched; the
      thin broken lines are the ends of stages 4 and 5.  Broken and dotted horizontal lines: the midpoints 6.7 and
      5.2 A that separate crystal-like from rotated-like frames.  A filled triangle marks the conversion by the
      rule fixed beforehand (first frame of the first visit of at least 1 ns in the crystal-like state).
(j, k) The two sensitivity runs of start R (zinc set 2s; candidate terms at N2).  The triangle of the second is
      open: its verdict is MIXED, and the triangle marks the start of its first crystal-like visit of at least
      1 ns, not a conversion.
Data: md/results/series/ (rule 6 of md/PLAN.md); verdicts, midpoints and times: facts md.* (section 40).
Drawing code: figS_md_common.py.
"""
import figS_md_common as C

if __name__ == '__main__':
    C.make('orientation', 'figS_md_orientation')
