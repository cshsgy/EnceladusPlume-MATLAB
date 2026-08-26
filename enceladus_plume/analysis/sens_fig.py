import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
d=np.load('/tmp/sensitivity.npz'); A,B,win,W=d['A'],d['B'],d['win'],d['winners']
# columns: 0 dw,1 L,2 phi2,3 alpha,4 sigma,5 chi2
mB=win==1; mA=win==0
fig,ax=plt.subplots(1,2,figsize=(9.2,3.8))
# (a) phi2 vs L, winning mode per realization
ax[0].scatter(W[mB,2],W[mB,1],c='tab:red',s=34,label='mode B wins (%d)'%mB.sum(),zorder=3,edgecolor='k',lw=0.3)
ax[0].scatter(W[mA,2],W[mA,1],c='tab:orange',s=34,marker='s',label='mode A wins (%d)'%mA.sum(),zorder=3,edgecolor='k',lw=0.3)
ax[0].scatter([123],[8.0],c='tab:red',marker='*',s=260,edgecolor='k',lw=0.6,zorder=4)
ax[0].scatter([84],[7.0],c='tab:orange',marker='*',s=260,edgecolor='k',lw=0.6,zorder=4)
ax[0].axhspan(4,12,color='tab:blue',alpha=0.08,zorder=0)
ax[0].text(140,11.4,'Hemingway & Mittal\n(2019): 4--12 km',fontsize=7,color='tab:blue',ha='center')
ax[0].set_xlabel(r'$2f$ phase $\phi_2$ [deg]'); ax[0].set_ylabel('depth $L$ [km]')
ax[0].set_xlim(20,170); ax[0].set_ylim(4,13); ax[0].legend(fontsize=8,loc='upper left')
ax[0].set_title('(a) depth stays shallow in both modes',fontsize=9)
ax[0].grid(alpha=0.3)
# (b) alpha histogram (winning mode)
ax[1].hist(W[:,3],bins=np.arange(0.2,1.45,0.08),color='0.6',edgecolor='k')
ax[1].axvline(0.81,color='tab:red',lw=2,label=r'adopted $\alpha=0.81$')
ax[1].set_xlim(0.2,1.45); ax[1].set_xlabel(r'$2f$ amplitude $\alpha$ (2f/1f)'); ax[1].set_ylabel('realizations')
ax[1].set_title(r'(b) strong $2f$ preferred (median $\alpha=0.82$)',fontsize=9)
ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
fig.tight_layout()
out='/home/sam2/dev/enceladus_plume_paper/Figures/sensitivity.pdf'
fig.savefig(out); print("wrote",out)
