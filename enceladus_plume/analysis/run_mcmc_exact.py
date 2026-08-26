"""Rigorous posterior: affine-invariant MCMC against the EXACT forward model
(no emulator), parallelized with a process pool (emcee-style split move)."""
import numpy as np, time
import multiprocessing as mp
import run_fit as R
from enceladus_plume.gas_dynamics.lookup import GasLookupTable

_ARGS=_LO=_HI=None
def _setup():
    global _ARGS,_LO,_HI
    lut=GasLookupTable('/tmp/encfig_lut.npz',clean=True); cfg=R._cfg()
    ma_o,y_o,sig=np.loadtxt(R._DATA,delimiter=',',skiprows=1).T
    weff_of=R.build_weff_interp(cfg)
    _ARGS=(weff_of,lut,cfg,ma_o,y_o,R._weights(ma_o,sig))
    _LO=np.array([b[0] for b in R.MLE_BOUNDS]); _HI=np.array([b[1] for b in R.MLE_BOUNDS])

def logprob(theta):
    th=np.asarray(theta,float)
    if np.any(th<_LO) or np.any(th>_HI): return -np.inf
    return -R._neg_loglike(th,*_ARGS)

def main():
    _setup()
    nd,nw,nsteps,a = 5,64,700,2.0
    rng=np.random.default_rng(0)
    mle=np.array([9.3,6.9,82.0,0.87,22.0])
    pos=np.clip(mle+rng.normal(0,[1.0,1.0,10.0,0.15,3.0],size=(nw,nd)),_LO+1e-6,_HI-1e-6)
    pool=mp.Pool(32)
    lp=np.array(pool.map(logprob,list(pos)))
    half=nw//2; idx=np.arange(nw); chain=np.zeros((nsteps,nw,nd)); t0=time.time()
    for s in range(nsteps):
        for k in (0,1):
            S=idx[k*half:(k+1)*half]; C=idx[(1-k)*half:(2-k)*half]
            props=[]; zs=np.empty(half)
            for m,i in enumerate(S):
                j=C[rng.integers(half)]; z=((a-1.0)*rng.random()+1.0)**2/a
                props.append(pos[j]+z*(pos[i]-pos[j])); zs[m]=z
            lpp=np.array(pool.map(logprob,props))
            acc=np.log(rng.random(half)) < (nd-1)*np.log(zs)+lpp-lp[S]
            for m,i in enumerate(S):
                if acc[m]: pos[i]=props[m]; lp[i]=lpp[m]
        chain[s]=pos
        if (s+1)%25==0:
            el=(time.time()-t0)/60
            np.save('/tmp/mcmc_exact_chain.npy',chain[:s+1])
            print("step %d/%d meanlogP=%.1f  %.1fmin  ETA %.0fmin"%(
                s+1,nsteps,lp.mean(),el,el/(s+1)*(nsteps-s-1)),flush=True)
    pool.close(); pool.join()
    samples=chain[nsteps//3:].reshape(-1,nd)
    q=np.percentile(samples,[16,50,84],axis=0)
    np.savez('/tmp/diurnal_posterior_exact.npz',samples=samples,med=q[1],lo=q[1]-q[0],up=q[2]-q[1],
             labels=["dw [mm]","L [km]","phi2 [deg]","alpha (2f/1f)","sigma [deg]"])
    print("EXACT MCMC DONE  median:",np.round(q[1],3),flush=True)

if __name__=="__main__":
    main()
