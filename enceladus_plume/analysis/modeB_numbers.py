import numpy as np, run_fit as R
from enceladus_plume.gas_dynamics.lookup import GasLookupTable
lut=GasLookupTable('/tmp/encfig_lut.npz',clean=True); cfg=R._cfg()
ma_o,y_o,sig=np.loadtxt(R._DATA,delimiter=',',skiprows=1).T
args=(R.build_weff_interp(cfg),lut,cfg,ma_o,y_o,R._weights(ma_o,sig)); dof=len(ma_o)-5
b=dict(np.load('/tmp/modeB.npz'))
dw,L,phi2,al,s=b['dw']*1e3,b['L']/1e3,b['harm_phase'],b['harm_scale'],b['sigma']
def chi2(sig_): return 2*R._neg_loglike([dw,L,phi2,al,sig_],*args)/dof
print("mode B ensemble (sigma=%.0f): chi2/dof=%.2f"%(s, chi2(s)))
print("mode B single   (sigma=0):  chi2/dof=%.2f"%(chi2(0.0)))
# exact posterior L interval for FITLPOST
p=np.load('/tmp/diurnal_posterior_exact.npz')['samples']
q=np.percentile(p[:,1],[16,50,84]); print("exact posterior L: %.1f (+%.1f/-%.1f)"%(q[1],q[2]-q[1],q[1]-q[0]))
