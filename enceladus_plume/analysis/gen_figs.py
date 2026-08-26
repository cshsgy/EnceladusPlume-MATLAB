import numpy as np, run_fit as R
from enceladus_plume.gas_dynamics.lookup import GasLookupTable
lut=GasLookupTable('/tmp/encfig_lut.npz',clean=True); cfg=R._cfg()
res=dict(np.load('/tmp/diurnal_fit_mle.npz'))   # mode B
# Fig 10 ensemble
R.plot_ensemble(res, lut, out='/tmp/diurnal_fit_ensemble.pdf')
# Fig 11 corner from EXACT posterior, truth-cross at mode B
P=np.load('/tmp/diurnal_posterior_exact.npz')
truths=[res['dw']*1e3, res['L']/1e3, res['harm_phase'], res['harm_scale'], res['sigma']]
R.plot_corner(P['samples'], list(P['labels']), truths=truths, out='/tmp/diurnal_fit_posterior.pdf')
print("truths(modeB)=", ["%.2f"%t for t in truths])
print("FIGS DONE")
