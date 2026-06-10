"""
make_figures.py
===============
Aggregate multi-seed results (experiments/results/MASTER_SUMMARY.csv + per-seed
CSVs) into publication figures with mean +/- std error bars over seeds.

Outputs -> experiments/results/figures/
"""
import os, glob
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(HERE, 'experiments', 'results')
FIG = os.path.join(RES, 'figures'); os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':7,
    'xtick.labelsize':6,'ytick.labelsize':6,'legend.fontsize':6,'legend.frameon':False,
    'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':400,'savefig.dpi':400,
    'savefig.bbox':'tight','pdf.fonttype':42})
PKD=r'pK$_d$'; C={'meanpool':'#95a5a6','allcdr':'#e74c3c','shuffle':'#bdc3c7'}
def save(fig,n):
    for e in ('pdf','png'): fig.savefig(os.path.join(FIG,f'{n}.{e}'),bbox_inches='tight')
    print('saved',n)
def panel(ax,l): ax.text(0.02,0.97,l,transform=ax.transAxes,fontsize=10,fontweight='bold',
    va='top',bbox=dict(boxstyle='round,pad=0.2',fc='white',alpha=0.7))

def agg(name):
    p=os.path.join(RES,name,'aggregated_summary.csv')
    return pd.read_csv(p) if os.path.isfile(p) else None

mp_cv   = agg('ours_meanpool_cv')
cdr_cv  = agg('ours_allcdr_cv')
mp_bn   = agg('ours_meanpool_benchmark')
cdr_bn  = agg('ours_allcdr_benchmark')
shuf_cv = agg('antigen_shuffle_cv')

# ── Fig A: mean-pool vs All-CDR CV across datasets (multi-seed) ──────────────
if mp_cv is not None and cdr_cv is not None:
    dss=['sabdab','abbind','skempi']; lbl=['SAbDab','AbBind','SKEMPI']
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.0))
    for ax,metric,ylab,letter in [(axes[0],'pearson','Pearson r','A'),(axes[1],'spearman','Spearman ρ','B')]:
        x=np.arange(3); w=0.38
        mv=[float(mp_cv[mp_cv.dataset==d][f'{metric}_mean'].values[0]) for d in dss]
        ms=[float(mp_cv[mp_cv.dataset==d][f'{metric}_std'].values[0]) for d in dss]
        cv_=[float(cdr_cv[cdr_cv.dataset==d][f'{metric}_mean'].values[0]) for d in dss]
        cs=[float(cdr_cv[cdr_cv.dataset==d][f'{metric}_std'].values[0]) for d in dss]
        ax.bar(x-w/2,mv,w,yerr=ms,color=C['meanpool'],alpha=0.85,label='Mean-pool',error_kw=dict(lw=1,capsize=3))
        ax.bar(x+w/2,cv_,w,yerr=cs,color=C['allcdr'],alpha=0.85,label='All-CDR',error_kw=dict(lw=1,capsize=3))
        ax.set_xticks(x); ax.set_xticklabels(lbl); ax.set_ylabel(ylab)
        ax.set_title(f'{letter}  {ylab} (mean ± std, 5 seeds)',fontsize=8); panel(ax,letter); ax.legend()
        for i in range(3):
            ax.text(i-w/2,mv[i]+ms[i]+0.01,f'{mv[i]:.2f}',ha='center',fontsize=5.5)
            ax.text(i+w/2,cv_[i]+cs[i]+0.01,f'{cv_[i]:.2f}',ha='center',fontsize=5.5)
    plt.tight_layout(w_pad=2); save(fig,'msfig1_cv_meanpool_vs_allcdr'); plt.close()

# ── Fig B: benchmark mean-pool vs All-CDR (multi-seed) ───────────────────────
if mp_bn is not None and cdr_bn is not None:
    fig,ax=plt.subplots(figsize=(3.6,3.2))
    methods=['Mean-pool','All-CDR']; cols=[C['meanpool'],C['allcdr']]
    vals=[float(mp_bn.pearson_mean.values[0]),float(cdr_bn.pearson_mean.values[0])]
    errs=[float(mp_bn.pearson_std.values[0]),float(cdr_bn.pearson_std.values[0])]
    ax.bar(range(2),vals,yerr=errs,color=cols,alpha=0.85,width=0.55,error_kw=dict(lw=1,capsize=4))
    ax.axhline(0.60,color='k',ls=':',lw=1,alpha=0.6,label='~0.60 reference')
    ax.set_xticks(range(2)); ax.set_xticklabels(methods); ax.set_ylabel('Benchmark Pearson r')
    ax.set_title('SAbDab→Benchmark (5 seeds)\nfull-epoch protocol',fontsize=8); ax.legend()
    for i,v in enumerate(vals): ax.text(i,v+errs[i]+0.01,f'{v:.3f}',ha='center',fontsize=7,fontweight='bold')
    plt.tight_layout(); save(fig,'msfig2_benchmark'); plt.close()

# ── Fig C: antigen importance (real vs shuffled, multi-seed) ─────────────────
if mp_cv is not None and shuf_cv is not None:
    dss=['sabdab','abbind','skempi']; lbl=['SAbDab','AbBind','SKEMPI']
    fig,ax=plt.subplots(figsize=(4.2,3.0)); x=np.arange(3); w=0.38
    rv=[float(mp_cv[mp_cv.dataset==d]['pearson_mean'].values[0]) for d in dss]
    re_=[float(mp_cv[mp_cv.dataset==d]['pearson_std'].values[0]) for d in dss]
    sv=[float(shuf_cv[shuf_cv.dataset==d]['pearson_mean'].values[0]) for d in dss]
    se=[float(shuf_cv[shuf_cv.dataset==d]['pearson_std'].values[0]) for d in dss]
    ax.bar(x-w/2,rv,w,yerr=re_,color=C['allcdr'],alpha=0.85,label='Real antigen',error_kw=dict(lw=1,capsize=3))
    ax.bar(x+w/2,sv,w,yerr=se,color=C['shuffle'],alpha=0.85,label='Shuffled antigen',error_kw=dict(lw=1,capsize=3))
    for i in range(3):
        d=rv[i]-sv[i]; ax.annotate(f'Δ{d:+.2f}',(i,max(rv[i],sv[i])+0.06),ha='center',fontsize=6.5,
            color=C['allcdr'] if d>0.05 else 'gray',fontweight='bold')
    ax.set_xticks(x); ax.set_xticklabels(lbl); ax.set_ylabel('Pearson r')
    ax.set_title('Antigen importance (5 seeds): real vs shuffled',fontsize=8); ax.legend()
    plt.tight_layout(); save(fig,'msfig3_antigen_shuffle'); plt.close()

# ── Fig D: master overview ───────────────────────────────────────────────────
master=os.path.join(RES,'MASTER_SUMMARY.csv')
if os.path.isfile(master):
    m=pd.read_csv(master)
    print('\nMASTER SUMMARY:'); print(m.to_string(index=False))
print('\nFigures ->',FIG)
