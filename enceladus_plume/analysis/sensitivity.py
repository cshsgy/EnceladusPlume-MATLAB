import numpy as np, run_fit as R, time
import multiprocessing as mp
from scipy.optimize import minimize
from enceladus_plume.gas_dynamics.lookup import GasLookupTable

SEEDS = {'A':[9.3,6.9,82.0,0.87,22.0], 'B':[6.77,7.83,122.2,0.81,32.17]}
_lut=_cfg=_weff=_ma=_sig=None

def _init():
    global _lut,_cfg,_weff,_ma,_sig
    _lut=GasLookupTable('/tmp/encfig_lut.npz',clean=True)
    _cfg=R._cfg(); _weff=R.build_weff_interp(_cfg, verbose=False)
    _ma,_,_sig=np.loadtxt(R._DATA,delimiter=',',skiprows=1).T

def _one(task):
    i, yp = task
    w=R._weights(_ma,_sig); args=(_weff,_lut,_cfg,_ma,yp,w)
    res={}
    for tag,x0 in SEEDS.items():
        loc=minimize(lambda th: R._neg_loglike(np.asarray(th),*args), x0,
                     method='Nelder-Mead', bounds=R.MLE_BOUNDS,
                     options=dict(xatol=1e-2,fatol=1e-2,maxiter=25))
        x=loc.x
        res[tag]=[float(x[0]),float(x[1]),float(x[2]),float(x[3]),float(x[4]),
                  float(2*R._neg_loglike(x,*args))]
    return i,res

if __name__=='__main__':
    ma,y0,sig=np.loadtxt(R._DATA,delimiter=',',skiprows=1).T
    np.random.seed(0); N=30
    tasks=[(i, np.clip(y0+np.random.normal(0,sig), 1.0, None)) for i in range(N)]
    t0=time.time()
    with mp.Pool(16, initializer=_init) as pool:
        out=pool.map(_one, tasks)
    A=np.array([r['A'] for _,r in out]); B=np.array([r['B'] for _,r in out])
    # winning mode per realization (lower chi2)
    win = np.where(B[:,5] <= A[:,5], 1, 0)  # 1 = mode B wins
    winners = np.where(win[:,None]==1, B, A)
    np.savez('/tmp/sensitivity.npz', A=A, B=B, win=win, winners=winners,
             labels=['dw','L','phi2','alpha','sigma','chi2'])
    lab=['dw[mm]','L[km]','phi2','alpha','sigma']
    def stat(arr,j): 
        q=np.percentile(arr[:,j],[16,50,84]); return q[1],q[2]-q[1],q[1]-q[0]
    print("N=%d realizations, %.1f min"%(N,(time.time()-t0)/60))
    print("mode B wins in %d/%d (%.0f%%); mode A wins in %d"%(win.sum(),N,100*win.mean(),N-win.sum()))
    print("--- WINNING-MODE parameter spread (median +/- 68%%) ---")
    for j,l in enumerate(lab): 
        m,up,lo=stat(winners,j); print("  %-7s %.2f +%.2f/-%.2f"%(l,m,up,lo))
    print("  chi2/dof %.2f +/- %.2f"%(np.median(winners[:,5])/(len(ma)-5), winners[:,5].std()/(len(ma)-5)))
    print("--- MODE B across all realizations ---")
    for j,l in enumerate(lab):
        m,up,lo=stat(B,j); print("  %-7s %.2f +%.2f/-%.2f"%(l,m,up,lo))
    print("--- MODE A across all realizations ---")
    for j,l in enumerate(lab):
        m,up,lo=stat(A,j); print("  %-7s %.2f +%.2f/-%.2f"%(l,m,up,lo))
    print("SENSITIVITY DONE", flush=True)
