import numpy as np, run_fit as R, time
from scipy.optimize import minimize
from enceladus_plume.gas_dynamics.lookup import GasLookupTable
lut=GasLookupTable('/tmp/encfig_lut.npz',clean=True); cfg=R._cfg()
ma_o,y_o,sig=np.loadtxt(R._DATA,delimiter=',',skiprows=1).T
w_o=R._weights(ma_o,sig); weff_of=R.build_weff_interp(cfg)
args=(weff_of,lut,cfg,ma_o,y_o,w_o); dof=len(ma_o)-5

def refine(x0, tag):
    t0=time.time()
    loc=minimize(lambda th: R._neg_loglike(np.asarray(th),*args), x0, method='Nelder-Mead',
                 bounds=R.MLE_BOUNDS, options=dict(xatol=1e-2,fatol=1e-2,maxiter=30))
    x=loc.x; chi2=2*R._neg_loglike(x,*args)
    dw,L,phi2,al,s=x; we=weff_of(dw,L); cfg.physical.equilibrium_depth=L*1e3
    MA,fl=R._flux_curve(cfg,L*1e3,dw*1e-3,we,lut,harm_scale=al,harm_phase=phi2)
    o=np.argsort(MA); g,fs=R._ensemble_smooth(MA[o],fl[o],s); p0,A,_=R._best_phi_A(g,fs,ma_o,y_o,w_o)
    res=dict(dw=dw*1e-3,L=L*1e3,w_eff=we,harm_scale=float(al),harm_phase=float(phi2),
             sigma=float(s),phi0=float(p0),A=float(A),chi2=float(chi2),dof=int(dof),chi2_red=float(chi2/dof))
    np.savez('/tmp/mode%s.npz'%tag,**res)
    print("mode %s [%.1fmin]: dw=%.1fmm L=%.1fkm phi2=%.0f alpha=%.2f sigma=%.0f  chi2/dof=%.2f"%(
        tag,(time.time()-t0)/60, dw,L,phi2,al,s, chi2/dof), flush=True)
    R.plot_overlay(res, lut, out='/tmp/mode%s_fit.pdf'%tag)
    return res

refine([9.3,6.9,82.0,0.87,22.0], "A")   # around the old MLE (phi2~82)
refine([6.77,7.83,122.2,0.81,32.17], "B")  # around the posterior mode (phi2~122)
print("BOTH MODES DONE", flush=True)
