# throwaway (version 1.3): a test protocol for two short runs of the CURRENT driver on its own stand-in systems:
# md/protocol.json as it stands (hydration check: six ring carbons, 0.25 nm) with purpose "test", output every 0.1 ps
# for the state log, the solute frames and the full frames, and at most 30 minimiser iterations per pass.  With
# --equil-scale 0.001 every step is then a whole number of 0.1 ps, so that each step ends on a full frame and the
# driver's count at the end of a step can be set beside the count of the analysis in that frame.
# usage: make_protocol.py MD_DIR OUT.json
import hashlib, json, sys
src = sys.argv[1] + '/protocol.json'
p = json.load(open(src))
p['purpose'] = 'test'
p['output']['state_log_interval_ps'] = 0.1
p['output']['solute_frame_interval_ps'] = 0.1
p['output']['full_frame_interval_ps'] = 0.1
p['minimisation']['max_iterations'] = 30
p['test_note'] = ('md/protocol.json (sha256 ' + hashlib.sha256(open(src, 'rb').read()).hexdigest() + ') with purpose "test", output '
                  'intervals of 0.1, 0.1 and 0.1 ps and at most 30 minimiser iterations per pass; nothing else is changed. '
                  'For the check of md/analysis/analyse_md.py (specification 1.3) against the records of the driver.')
json.dump(p, open(sys.argv[2], 'w'), indent=1)
print(json.dumps(p['equilibration']['hydration_check']['site_atoms']), p['equilibration']['hydration_check']['water_cutoff_nm'])
