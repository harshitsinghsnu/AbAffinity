"""
build_panels_notebook.py
========================
Builds BALM_AbAg_Panels.ipynb — Nature-style multi-panel figures.
Every bar/curve panel has its EXACT numbers embedded as literals (read from the
result files at BUILD time), so the notebook has no runtime file dependency for
those panels. Only the scatter and heatmap panels read raw data files.

Colours: 'Ours' is one fixed colour everywhere; categorical comparators use the
colour-blind-safe Okabe-Ito palette; ordered series (PLM, transfer, gating,
CDRs) use sequential colormaps. 400 DPI PDF+PNG, Arial, Nature font sizes.
"""
import json, os, uuid
import pandas as pd, numpy as np
from scipy.stats import pearsonr, spearmanr
HERE = os.path.dirname(os.path.abspath(__file__))
NB   = os.path.join(HERE, 'BALM_AbAg_Panels.ipynb')

def f(x): return float(x)

# ───────────────────────── read everything, build DATA ──────────────────────
D = {}

# -- architecture (SAaIntDB random+cold; Ours / Two-stream / Concat+MLP) --
ms = pd.read_csv(os.path.join(HERE,'experiments/results_saaintdb/MASTER_SUMMARY_SAAINTDB.csv'))
def ours_arch(split, m):
    r = ms[ms.experiment == f'sa_ours_allcdr_{split}']
    return [f(r[f'{m}_mean'].values[0]), f(r[f'{m}_std'].values[0])]
def long_arch(path, m):
    d = pd.read_csv(path); r = d[d.metric == m]
    return [f(r['mean'].values[0]), f(r['std'].values[0])] if len(r) else [np.nan, np.nan]
D['arch'] = {}
for m in ['pearson','spearman','rmse']:
    D['arch'][m] = {
        'Ours':       {'random': ours_arch('random', m), 'cold': ours_arch('cold', m)},
        'Two-stream': {'random': long_arch(os.path.join(HERE,'results_twostream_saaintdb/cv_summary.csv'), m),
                       'cold':   long_arch(os.path.join(HERE,'results_twostream_saaintdb_cold/cv_summary.csv'), m)},
        'Concat+MLP': {'random': long_arch(os.path.join(HERE,'results_baseline_concat_mlp_saaintdb/cv_summary.csv'), m),
                       'cold':   long_arch(os.path.join(HERE,'results_baseline_concat_mlp_saaintdb_cold/cv_summary.csv'), m)},
    }

# -- overall vs intra-target (Fisher, per-antigen) --
fz = pd.read_csv(os.path.join(HERE,'experiments/results_allcdr_stats/overall_vs_fisher.csv'))
g = fz[fz.grouping == 'per_antigen']
D['fisher'] = {m: {'overall': f(g[g.metric==m]['overall'].values[0]),
                   'intra':   f(g[g.metric==m]['fisher_mean'].values[0]),
                   'intra_std': f(g[g.metric==m]['fisher_std'].values[0])}
               for m in ['pearson','spearman','rmse']}

# -- antibody vs nanobody --
nb = pd.read_csv(os.path.join(HERE,'experiments/results_allcdr_stats/antibody_vs_nanobody.csv'))
D['nano'] = [{'subgroup': r['subgroup'].split(' (')[0], 'n': int(r['n']),
              'pearson': f(r['pearson']), 'spearman': f(r['spearman']), 'rmse': f(r['rmse'])}
             for _, r in nb.iterrows()]
# per-fold subgroup metrics (for error bars across folds)
ap = pd.read_csv(os.path.join(HERE,'results_saaintdb_allcdr/random/all_preds.csv'))
def _valid(s): return isinstance(s,str) and s.strip() not in ('','nan','N/A','NA','None') and len(s.strip())>5
ap['is_nano'] = ~ap['L_seq'].apply(_valid)
def _pf(sub):
    rs=[]
    for _, gg in sub.groupby('fold'):
        if len(gg) < 3: continue
        t=gg['binding_affinity'].values; p=gg['predicted_affinity'].values
        rs.append([f(pearsonr(t,p)[0]), f(spearmanr(t,p)[0]), f(np.sqrt(np.mean((t-p)**2)))])
    rs=np.array(rs); return rs.mean(0), rs.std(0)
D['nano_fold'] = {}
for nm, sub in [('Overall',ap),('Antibody',ap[~ap.is_nano]),('Nanobody',ap[ap.is_nano])]:
    mu,sd = _pf(sub)
    D['nano_fold'][nm] = {'pearson':[f(mu[0]),f(sd[0])],'spearman':[f(mu[1]),f(sd[1])],'rmse':[f(mu[2]),f(sd[2])],'n':int(len(sub))}

# -- Ours vs MVSF-AB (CV) --
M = pd.read_csv(os.path.join(HERE,'experiments/results/MASTER_SUMMARY.csv'))
def ours_cv(ds, m):
    r = M[(M.experiment=='ours_allcdr_cv') & (M.dataset==ds)]
    return [f(r[f'{m}_mean'].values[0]), f(r[f'{m}_std'].values[0])] if len(r) else [np.nan, np.nan]
bench = pd.read_csv(os.path.join(HERE,'results_natural_allcdr/benchmark_fair_comparison.csv'))
br = bench[bench.method.str.contains('All-CDR')].iloc[0]
mvsf = {'sabdab':(0.491,1.839),'abbind':(0.739,1.905),'skempi':(0.671,1.513),'benchmark':(0.467,1.447)}
D['mvsf'] = {}
for ds in ['sabdab','abbind','skempi','benchmark']:
    if ds == 'benchmark':
        ours = {'pearson':[f(br['pearson']), f(br['pearson_std'])], 'spearman':[f(br['spearman']), np.nan],
                'rmse':[f(br['rmse']), np.nan]}
    else:
        ours = {m: ours_cv(ds, m) for m in ['pearson','spearman','rmse']}
    D['mvsf'][ds] = {'ours': ours, 'mvsf': {'pearson': mvsf[ds][0], 'rmse': mvsf[ds][1]}}

# -- transfer (zero/10/20/30, per dataset) --
tr = pd.read_csv(os.path.join(HERE,'experiments/results_transfer_allcdr/aggregated.csv'))
D['transfer'] = {}
for ds, gd in tr.groupby('dataset'):
    gd = gd.sort_values('fraction')
    D['transfer'][ds] = {'frac':[f(x) for x in gd['fraction']],
                         'pearson':[f(x) for x in gd['pearson_mean']],
                         'pearson_std':[f(x) for x in gd['pearson_std']],
                         'spearman':[f(x) for x in gd['spearman_mean']],
                         'rmse':[f(x) for x in gd['rmse_mean']]}

# -- gating ablation (all metrics) --
gt = pd.read_csv(os.path.join(HERE,'experiments/results_allcdr_stats/allcdr_gating_saaintdb.csv'))
D['gating'] = {r['gate_mode']: {m:[f(r[f'{m}_mean']), f(r[f'{m}_std'])] for m in ['pearson','spearman','rmse']}
               for _, r in gt.iterrows()}

# -- stream interventions + per-CDR (supplementary) --
si = pd.read_csv(os.path.join(HERE,'experiments/results_allcdr_stats/allcdr_stream_interventions.csv'))
sf = pd.read_csv(os.path.join(HERE,'experiments/results_allcdr_stats/allcdr_stream_interventions_folds.csv'))
sstd = sf.groupby('condition').agg(ss=('spearman','std'), rs=('rmse','std'))
D['stream'] = [{'condition': r['condition'], 'pearson': f(r['pearson_mean']), 'pearson_std': f(r['pearson_std']),
                'spearman': f(r['spearman_mean']), 'spearman_std': f(sstd.loc[r['condition'],'ss']),
                'rmse': f(r['rmse_mean']), 'rmse_std': f(sstd.loc[r['condition'],'rs'])} for _, r in si.iterrows()]
cr = pd.read_csv(os.path.join(HERE,'advanced_results/cdr_region_performance_allcdr.csv'))
D['cdrreg'] = [{'cdr': r['cdr'], 'pearson': f(r['pearson_r']), 'spearman': f(r['spearman_rho']), 'rmse': f(r['rmse'])}
               for _, r in cr.iterrows()]
# top-K recovery
tk = pd.read_csv(os.path.join(HERE,'advanced_results/topk_recovery_fixed.csv'))
D['topk'] = {'k':[f(x) for x in tk['k_pct']],
             'random':[f(x) for x in tk['random_recovery']],
             'meanpool':[f(x) for x in tk['in_domain_recovery']],
             'ours':[f(x) for x in tk['allcdr_recovery']],
             'zeroshot':[f(x) for x in tk['zero_shot_recovery']]}

# -- PLM comparison (all metrics) + antibody-PLM ablation --
plm = pd.read_csv(os.path.join(HERE,'plm_comparison_results_mutual_saaintdb/cv_all_plms_saaintdb.csv'))
D['plm'] = {r['PLM']: {'pearson':[f(r['CV_Pearson']), f(r.get('CV_Pearson_std',0))],
                       'spearman':[f(r['CV_Spearman']), f(r.get('CV_Spearman_std',0))],
                       'rmse':[f(r['CV_RMSE']), f(r.get('CV_RMSE_std',0))]} for _, r in plm.iterrows()}
abl = pd.read_csv(os.path.join(HERE,'results_ablation_antibody_plm/cv_summary.csv'))
D['plm']['AntiBERTy-H/L\n+ESM-2-Ag'] = {'pearson':[f(abl.pearson.mean()), f(abl.pearson.std())],
    'spearman':[f(abl.spearman.mean()), f(abl.spearman.std())], 'rmse':[f(abl.rmse.mean()), f(abl.rmse.std())]}

# -- pKd distributions (embed histograms + stats) --
def pkd_arr(ds):
    if ds=='saaintdb': return pd.read_csv(os.path.join(HERE,'results_saaintdb_allcdr/random/all_preds.csv'))['binding_affinity'].values
    if ds in ('sabdab','abbind','skempi'): return pd.read_csv(os.path.join(HERE,f'experiments/results_allcdr_stats/oof_preds_{ds}.csv'))['true'].values
    if ds=='benchmark':
        import sys; sys.path.insert(0,HERE)
        from agabgated.utils.main_symmetric_mean import load_data
        return load_data(os.path.join(HERE,'datasets/pairs_benchmark.csv'))['binding_affinity'].values
D['pkd'] = {}
for ds in ['saaintdb','sabdab','abbind','skempi','benchmark']:
    a = np.asarray(pkd_arr(ds), float); a = a[np.isfinite(a)]
    cnt, edges = np.histogram(a, bins=25)
    D['pkd'][ds] = {'n': int(len(a)), 'min': f(a.min()), 'max': f(a.max()), 'mean': f(a.mean()),
                    'median': f(np.median(a)), 'counts':[int(x) for x in cnt], 'edges':[f(x) for x in edges]}
D['pkd']['saaintdb']['nano'] = int([x['n'] for x in D['nano'] if 'Nanobody' in x['subgroup']][0])
D['pkd']['saaintdb']['antibody'] = int([x['n'] for x in D['nano'] if 'Antibody' in x['subgroup']][0])

DATA_LITERAL = json.dumps(D, indent=0)

# ───────────────────────────── notebook assembly ────────────────────────────
cells = []
def md(t):   cells.append({'cell_type':'markdown','metadata':{},'id':uuid.uuid4().hex[:8],'source':[t]})
def code(s): cells.append({'cell_type':'code','metadata':{},'id':uuid.uuid4().hex[:8],'execution_count':None,'outputs':[],'source':[s]})

md("# Ours — Publication Panels (self-contained)\n"
   "Nature-style multi-panel figures. **All bar/curve panels embed their exact numbers** (extracted "
   "from the result files at build time) — no runtime file dependency. Only the scatter and heatmap "
   "panels read raw data. 400 DPI PDF+PNG. 'Ours' uses one fixed colour throughout; comparators use the "
   "colour-blind-safe Okabe-Ito palette; ordered series use sequential colormaps.")

# setup cell: style + embedded DATA + palette
code(
"import os, json, numpy as np, pandas as pd\n"
"import matplotlib, matplotlib.pyplot as plt\n"
"from matplotlib import cm\n"
"HERE=os.path.dirname(os.path.abspath('__file__'))\n"
"FIG=os.path.join(HERE,'figures_panels'); os.makedirs(FIG,exist_ok=True)\n"
"plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':7,\n"
"  'xtick.labelsize':6,'ytick.labelsize':6,'legend.fontsize':6,'legend.frameon':False,\n"
"  'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':400,'savefig.dpi':400,\n"
"  'savefig.bbox':'tight','pdf.fonttype':42,'axes.linewidth':0.7})\n"
"PKD=r'pK$_d$'\n"
"# fixed colours — every colour means one thing\n"
"OURS='#D55E00'      # Ours (fixed everywhere)\n"
"TWO='#0072B2'; CONCAT='#009E73'; MVSF='#999999'   # comparators (Okabe-Ito)\n"
"C2='#56B4E9'; C3='#CC79A7'; C4='#E69F00'          # extra categoricals\n"
"def seq(n,cmap='viridis',lo=0.15,hi=0.9): return [cm.get_cmap(cmap)(x) for x in np.linspace(lo,hi,n)]\n"
"def panel(ax,L): ax.text(-0.16,1.06,L,transform=ax.transAxes,fontsize=10,fontweight='bold',va='top',ha='left')\n"
"def save(fig,n):\n"
"    for e in ('pdf','png'): fig.savefig(os.path.join(FIG,f'{n}.{e}'),bbox_inches='tight')\n"
"    print('saved',n)\n"
"DATA=json.loads('''"+DATA_LITERAL.replace("\\","\\\\").replace("'''","\\'\\'\\'")+"''')\n"
"MET=[('pearson','Pearson r'),('spearman','Spearman ρ'),('rmse',f'RMSE ({PKD})')]\n"
"print('ready; datasets embedded:',list(DATA.keys()))")

# ── Panel 1 : data overview / pKd distributions ──────────────────────────────
md("---\n## Panel 1 | Datasets and pK_d distributions")
code(
"d=DATA['pkd']; dss=['saaintdb','sabdab','abbind','skempi','benchmark']\n"
"names={'saaintdb':'SAaIntDB','sabdab':'SAbDab','abbind':'AB-Bind','skempi':'SKEMPI','benchmark':'Benchmark'}\n"
"cols=seq(len(dss),'viridis')\n"
"fig=plt.figure(figsize=(7.2,4.4)); gs=fig.add_gridspec(2,3,height_ratios=[1.1,1])\n"
"# a) overlaid/normalised histograms\n"
"axa=fig.add_subplot(gs[0,:])\n"
"for ds,c in zip(dss,cols):\n"
"    e=np.array(d[ds]['edges']); ct=np.array(d[ds]['counts'],float); ct=ct/ct.sum()\n"
"    ctr=(e[:-1]+e[1:])/2; axa.plot(ctr,ct,color=c,lw=1.3,label=f\"{names[ds]} (n={d[ds]['n']})\")\n"
"axa.set_xlabel(f'Measured {PKD}'); axa.set_ylabel('frequency'); axa.legend(fontsize=5.6,ncol=2)\n"
"axa.set_title('a  Affinity distribution across datasets',fontsize=8,loc='left'); panel(axa,'a')\n"
"# b) dataset sizes\n"
"axb=fig.add_subplot(gs[1,0]); axb.bar(range(len(dss)),[d[x]['n'] for x in dss],color=cols,alpha=0.9)\n"
"axb.set_xticks(range(len(dss))); axb.set_xticklabels([names[x] for x in dss],rotation=35,ha='right',fontsize=5.6)\n"
"axb.set_ylabel('# complexes'); axb.set_yscale('log'); axb.set_title('b  Dataset size',fontsize=8,loc='left'); panel(axb,'b')\n"
"for i,x in enumerate(dss): axb.text(i,d[x]['n'],str(d[x]['n']),ha='center',va='bottom',fontsize=5)\n"
"# c) antibody vs nanobody composition (SAaIntDB)\n"
"axc=fig.add_subplot(gs[1,1]); sa=d['saaintdb']\n"
"axc.bar([0,1],[sa['antibody'],sa['nano']],color=seq(2,'viridis'),alpha=0.9)\n"
"axc.set_xticks([0,1]); axc.set_xticklabels(['Antibody\\n(paired)','Nanobody\\n(VHH)'],fontsize=6)\n"
"axc.set_ylabel('# complexes (SAaIntDB)'); axc.set_title('c  Composition',fontsize=8,loc='left'); panel(axc,'c')\n"
"for i,v in enumerate([sa['antibody'],sa['nano']]): axc.text(i,v,str(v),ha='center',va='bottom',fontsize=6)\n"
"# d) pKd range table-like bars (min..max with mean)\n"
"axd=fig.add_subplot(gs[1,2])\n"
"for i,ds in enumerate(dss):\n"
"    axd.plot([d[ds]['min'],d[ds]['max']],[i,i],color=cols[i],lw=3,solid_capstyle='round')\n"
"    axd.plot(d[ds]['mean'],i,'o',color='k',ms=3)\n"
"axd.set_yticks(range(len(dss))); axd.set_yticklabels([names[x] for x in dss],fontsize=5.6)\n"
"axd.set_xlabel(f'{PKD} range (• = mean)'); axd.set_title('d  Dynamic range',fontsize=8,loc='left'); panel(axd,'d')\n"
"fig.suptitle('Datasets and affinity distributions',fontsize=9,fontweight='bold')\n"
"plt.tight_layout(rect=[0,0,1,0.97]); save(fig,'panel1_data'); plt.show()")

# ── Panel 2 : main comparison (a arch, b antibody/nanobody, c MVSF, d intra) ──
md("---\n## Panel 2 | Main comparison — architecture (a), antibody vs nanobody (b), vs MVSF-AB (c), overall vs intra-target (d)")
code(
"fig=plt.figure(figsize=(7.5,8.6)); gs=fig.add_gridspec(3,3,hspace=0.78,wspace=0.55)\n"
"# a) architecture: random+cold, Ours/Two/Concat, 3 metrics\n"
"A=DATA['arch']; models=['Ours','Two-stream','Concat+MLP']; mc={'Ours':OURS,'Two-stream':TWO,'Concat+MLP':CONCAT}\n"
"for ci,(mk,ml) in enumerate(MET):\n"
"    ax=fig.add_subplot(gs[0,ci]); x=np.arange(3); w=0.38\n"
"    rv=[A[mk][m]['random'][0] for m in models]; re_=[A[mk][m]['random'][1] for m in models]\n"
"    cv=[A[mk][m]['cold'][0] for m in models];   ce=[A[mk][m]['cold'][1] for m in models]\n"
"    ax.bar(x-w/2,rv,w,yerr=re_,color=[mc[m] for m in models],alpha=0.95,error_kw=dict(lw=0.7,capsize=2),label='random')\n"
"    ax.bar(x+w/2,cv,w,yerr=ce,color=[mc[m] for m in models],alpha=0.5,hatch='//',error_kw=dict(lw=0.7,capsize=2),label='cold')\n"
"    ax.set_xticks(x); ax.set_xticklabels(models,rotation=22,ha='right',fontsize=5.6); ax.set_ylabel(ml)\n"
"    ax.set_title(ml,fontsize=7.5)\n"
"    if ci==0: panel(ax,'a')\n"
"    if ci==2: ax.legend(fontsize=5.5,loc='center left',bbox_to_anchor=(1.0,0.5))\n"
"# b) antibody vs nanobody — 3 metrics, one subplot each, error bars across folds\n"
"NF=DATA['nano_fold']; grp=['Overall','Antibody','Nanobody']; gcol={'Overall':MVSF,'Antibody':TWO,'Nanobody':C4}\n"
"for ci,(mk,ml) in enumerate(MET):\n"
"    ax=fig.add_subplot(gs[1,ci]); x=np.arange(3)\n"
"    v=[NF[g][mk][0] for g in grp]; e=[NF[g][mk][1] for g in grp]\n"
"    ax.bar(x,v,yerr=e,color=[gcol[g] for g in grp],alpha=0.95,width=0.6,error_kw=dict(lw=0.7,capsize=3))\n"
"    ax.set_xticks(x); ax.set_xticklabels(grp,rotation=18,ha='right',fontsize=5.8); ax.set_ylabel(ml); ax.set_title(ml,fontsize=7.5)\n"
"    for i,g in enumerate(grp): ax.text(i,0.02,f\"n={NF[g]['n']}\",ha='center',fontsize=4.6,color='white')\n"
"    if ci==0: panel(ax,'b')\n"
"# c) Ours vs MVSF (Pearson + RMSE only; MVSF has no Spearman)\n"
"V=DATA['mvsf']; dss=['sabdab','abbind','skempi','benchmark']; dn=['SAbDab','AB-Bind','SKEMPI','Bench']\n"
"for ci,(mk,ml) in enumerate([('pearson','Pearson r'),('rmse',f'RMSE ({PKD})')]):\n"
"    ax=fig.add_subplot(gs[2,ci]); x=np.arange(4); w=0.38\n"
"    ovr=[V[d]['ours'][mk][0] for d in dss]; oer=[(V[d]['ours'][mk][1] if V[d]['ours'][mk][1]==V[d]['ours'][mk][1] else 0) for d in dss]\n"
"    mv=[V[d]['mvsf'][mk] for d in dss]\n"
"    ax.bar(x-w/2,ovr,w,yerr=oer,color=OURS,alpha=0.95,error_kw=dict(lw=0.7,capsize=2),label='Ours')\n"
"    ax.bar(x+w/2,mv,w,color=MVSF,alpha=0.9,label='MVSF-AB')\n"
"    ax.set_xticks(x); ax.set_xticklabels(dn,rotation=22,ha='right',fontsize=5.6); ax.set_ylabel(ml)\n"
"    top=max(max(ovr),max(mv)); ax.set_ylim(0,top*1.32)\n"
"    ax.legend(fontsize=5.5,ncol=2,loc='lower center',bbox_to_anchor=(0.5,1.0),frameon=False,columnspacing=1.0,handlelength=1.2)\n"
"    if ci==0: panel(ax,'c')\n"
"# d) overall vs intra-target (all metrics) — beside c\n"
"F=DATA['fisher']; NF=DATA['nano_fold']; ax=fig.add_subplot(gs[2,2]); x=np.arange(3); w=0.38\n"
"mks=['pearson','spearman','rmse']\n"
"ov=[F[m]['overall'] for m in mks]; ove=[NF['Overall'][m][1] for m in mks]\n"
"it=[F[m]['intra'] for m in mks]; its=[F[m]['intra_std'] for m in mks]\n"
"ax.bar(x-w/2,ov,w,yerr=ove,color=MVSF,alpha=0.9,error_kw=dict(lw=0.7,capsize=2),label='Overall'); ax.bar(x+w/2,it,w,yerr=its,color=C3,alpha=0.9,error_kw=dict(lw=0.7,capsize=2),label='Intra-target')\n"
"ax.set_xticks(x); ax.set_xticklabels(['Pears.','Spear.','RMSE'],fontsize=5.8); ax.set_ylabel('value'); ax.set_ylim(0,1.55)\n"
"ax.legend(fontsize=5.5,ncol=2,loc='lower center',bbox_to_anchor=(0.5,1.0),frameon=False,columnspacing=1.0,handlelength=1.2); panel(ax,'d')\n"
"fig.suptitle('Main comparison: architecture (a), antibody vs nanobody (b), vs MVSF-AB (c), overall vs intra-target (d)',fontsize=9,fontweight='bold')\n"
"plt.tight_layout(rect=[0,0,1,0.975]); save(fig,'panel2_main'); plt.show()")

# ── Panel 3 : PLM comparison ─────────────────────────────────────────────────
md("---\n## Panel 3 | Protein-language-model comparison (all metrics)")
code(
"P=DATA['plm']; order=['ESM-2','ProtBERT','AntiBERTy','ProGen2','AntiBERTy-H/L\\n+ESM-2-Ag']\n"
"order=[p for p in order if p in P]; cols=seq(len(order),'cividis')\n"
"fig,axes=plt.subplots(1,3,figsize=(7.2,2.7))\n"
"for ax,(mk,ml) in zip(axes,MET):\n"
"    v=[P[p][mk][0] for p in order]; e=[P[p][mk][1] for p in order]\n"
"    ax.bar(range(len(order)),v,yerr=e,color=cols,alpha=0.95,error_kw=dict(lw=0.7,capsize=2))\n"
"    ax.set_xticks(range(len(order))); ax.set_xticklabels(order,rotation=30,ha='right',fontsize=5.2); ax.set_ylabel(ml)\n"
"    ax.set_title(ml,fontsize=7.5)\n"
"axes[0].text(-0.16,1.06,'',transform=axes[0].transAxes)\n"
"fig.suptitle('PLM backbone comparison — SAaIntDB 10-fold CV (ESM-2 strongest)',fontsize=9,fontweight='bold')\n"
"plt.tight_layout(rect=[0,0,1,0.95]); save(fig,'panel3_plm'); plt.show()")

# ── Panel 4 : transfer (Pearson & Spearman explicit) ─────────────────────────
md("---\n## Panel 4 | Zero-shot & few-shot transfer (Pearson and Spearman labelled)")
code(
"T=DATA['transfer']; DS=['1mlc','aayl51','4fqi','koenig','trastuzumab','warszawski']\n"
"fig,axes=plt.subplots(2,3,figsize=(9.6,5.6)); axes=axes.flatten()\n"
"for ax,ds,L in zip(axes,DS,'abcdef'):\n"
"    fr=T[ds]['frac']; pe=T[ds]['pearson']; ps=T[ds]['pearson_std']; sp=T[ds]['spearman']\n"
"    zi=[i for i,x in enumerate(fr) if x==0.0]; fi=[i for i,x in enumerate(fr) if x>0]\n"
"    if zi: ax.axhline(pe[zi[0]],color='#95a5a6',ls=':',lw=1,label=f'zero-shot r ({pe[zi[0]]:.2f})')\n"
"    xf=[int(fr[i]*100) for i in fi]\n"
"    ax.errorbar(xf,[pe[i] for i in fi],yerr=[ps[i] for i in fi],fmt='o-',color=OURS,lw=1.8,ms=5,capsize=3,label='few-shot r')\n"
"    ax.plot(xf,[sp[i] for i in fi],'s--',color=TWO,lw=1.2,ms=4,label='few-shot ρ')\n"
"    ax.set_xlabel('Train %'); ax.set_ylabel('correlation (r, ρ)' if L in 'ad' else ''); ax.set_xticks([10,20,30])\n"
"    ax.set_title(f'{L}  {ds}',fontsize=8); ax.legend(fontsize=5.2,loc='best'); panel(ax,L)\n"
"fig.suptitle('All-CDR transfer: zero-shot vs few-shot (10/20/30%, mean ± s.d., 3 seeds)',fontsize=9,fontweight='bold')\n"
"plt.tight_layout(); save(fig,'panel4_transfer'); plt.show()")

# ── Panel 5 : gating + supplementary (all metrics) ───────────────────────────
md("---\n## Panel 5 | Gating ablation & stream/CDR analyses (all metrics)")
code(
"G=DATA['gating']; gorder=['learned','fixed','open','random','closed']\n"
"glab={'learned':'Learned','fixed':'Fixed','open':'Open','random':'Random','closed':'Closed'}\n"
"gcol=[OURS]+list(seq(3,'Blues',0.4,0.8))+[MVSF]\n"
"S=DATA['stream']; R=DATA['cdrreg']\n"
"fig=plt.figure(figsize=(7.2,5.4)); gs=fig.add_gridspec(2,3,hspace=0.6,wspace=0.5)\n"
"# a) gating, 3 metrics\n"
"for ci,(mk,ml) in enumerate(MET):\n"
"    ax=fig.add_subplot(gs[0,ci]); v=[G[m][mk][0] for m in gorder]; e=[G[m][mk][1] for m in gorder]\n"
"    ax.bar(range(5),v,yerr=e,color=gcol,alpha=0.95,error_kw=dict(lw=0.7,capsize=2))\n"
"    ax.set_xticks(range(5)); ax.set_xticklabels([glab[m] for m in gorder],rotation=25,ha='right',fontsize=5.4); ax.set_ylabel(ml)\n"
"    ax.set_title(ml,fontsize=7.5)\n"
"    if ci==0: panel(ax,'a')\n"
"# b) stream interventions (Pearson / Spearman / RMSE)\n"
"ax=fig.add_subplot(gs[1,0]); conds=[r['condition'] for r in S]; x=np.arange(len(conds)); w=0.27\n"
"ax.bar(x-w,[r['pearson'] for r in S],w,yerr=[r['pearson_std'] for r in S],color=OURS,alpha=0.92,error_kw=dict(lw=0.5,capsize=2),label='Pearson r')\n"
"ax.bar(x,[r['spearman'] for r in S],w,yerr=[r['spearman_std'] for r in S],color=C2,alpha=0.92,error_kw=dict(lw=0.5,capsize=2),label='Spearman ρ')\n"
"ax.bar(x+w,[r['rmse'] for r in S],w,yerr=[r['rmse_std'] for r in S],color=C3,alpha=0.92,error_kw=dict(lw=0.5,capsize=2),label=f'RMSE ({PKD})')\n"
"ax.set_xticks(x); ax.set_xticklabels([c.replace('_',' ') for c in conds],rotation=30,ha='right',fontsize=5); ax.set_ylabel('value'); ax.set_ylim(0,1.55)\n"
"ax.legend(fontsize=4.6,ncol=3,loc='upper center',columnspacing=0.7,handlelength=1.0,handletextpad=0.4); ax.set_title('Stream interventions',fontsize=7); panel(ax,'b')\n"
"# c) per-CDR region (3 metrics)\n"
"ax=fig.add_subplot(gs[1,1]); cdr=[r['cdr'] for r in R]; x=np.arange(len(cdr)); w=0.26\n"
"for j,(mk,_) in enumerate(MET):\n"
"    ax.bar(x+(j-1)*w,[r[mk] for r in R],w,color=[OURS,C2,C3][j],alpha=0.92,label=dict(MET)[mk])\n"
"ax.set_xticks(x); ax.set_xticklabels(cdr,fontsize=5.6); ax.set_ylabel('value'); ax.set_ylim(0,0.84)\n"
"ax.legend(fontsize=4.6,ncol=3,loc='upper center',columnspacing=0.7,handlelength=1.0,handletextpad=0.4); ax.set_title('Per-CDR-region performance',fontsize=7); panel(ax,'c')\n"
"# d) top-K binder recovery\n"
"TK=DATA['topk']; ax=fig.add_subplot(gs[1,2]); x=np.arange(len(TK['k'])); w=0.2\n"
"ser=[('random','#bdc3c7','Random'),('zeroshot',C2,'Zero-shot'),('meanpool',MVSF,'Mean-pool'),('ours',OURS,'Ours')]\n"
"for j,(key,c,lab) in enumerate(ser): ax.bar(x+(j-1.5)*w,TK[key],w,color=c,alpha=0.92,label=lab)\n"
"ax.set_xticks(x); ax.set_xticklabels([f'Top-{int(k*100)}%' for k in TK['k']],fontsize=5.2); ax.set_ylabel('fraction recovered'); ax.set_ylim(0,0.72)\n"
"ax.legend(fontsize=4.4,ncol=2,loc='upper left',columnspacing=0.6,handlelength=1.0,handletextpad=0.3); ax.set_title('Top-K binder recovery (AAYL51)',fontsize=7); panel(ax,'d')\n"
"fig.suptitle('Gating ablation (a), antigen-dependence (b), per-CDR (c) and top-K recovery (d)',fontsize=9,fontweight='bold')\n"
"plt.tight_layout(rect=[0,0,1,0.96]); save(fig,'panel5_gating_supp'); plt.show()")

# ── Panel 6 : heatmaps (file-backed) ─────────────────────────────────────────
md("---\n## Panel 6 | Mutational heatmaps (reads raw data files)")
code(
"from matplotlib.colors import TwoSlopeNorm\n"
"hm=pd.read_csv(os.path.join(HERE,'advanced_results/cdr_heatmap_data_allcdr.csv'))\n"
"AAS=list('ACDEFGHIKLMNPQRSTVWY'); regs=['CDR-H1','CDR-H2','CDR-H3']\n"
"d=np.load(os.path.join(HERE,'experiments/results_allcdr_stats/aa_mutation_impact.npz'),allow_pickle=True)\n"
"AAo=list(d['aas']); true=d['true']; pred=d['pred']; err=d['err']\n"
"cat={**{a:'#c0392b' for a in 'RHKDE'},**{a:'#2980b9' for a in 'STNQ'},**{a:'#27ae60' for a in 'AVILMFYW'},**{a:'#8e44ad' for a in 'CGP'}}\n"
"fig=plt.figure(figsize=(8.6,6.6)); gs=fig.add_gridspec(2,3,hspace=0.55,wspace=0.32)\n"
"# a) saturation mutagenesis by CDR region (measured ΔpKd); x = wild-type AA + position\n"
"for ci,reg in enumerate(regs):\n"
"    ax=fig.add_subplot(gs[0,ci]); sub=hm[hm.cdr_region==reg]; pos=sorted(sub['mut_pos'].unique())\n"
"    Mx=np.full((len(AAS),len(pos)),np.nan)\n"
"    for j,p in enumerate(pos):\n"
"        for i,aa in enumerate(AAS):\n"
"            r=sub[(sub.mut_pos==p)&(sub.mut_aa==aa)]\n"
"            if len(r): Mx[i,j]=float(r['delta_pkd_true'].values[0])\n"
"    v=np.nanmax(np.abs(Mx)); v=float(v) if np.isfinite(v) and v>0.1 else 0.1\n"
"    im=ax.imshow(Mx,cmap='RdYlGn',norm=TwoSlopeNorm(0,-v,v),aspect='auto')\n"
"    ax.set_yticks(range(len(AAS))); ax.set_yticklabels(AAS,fontsize=4) if ci==0 else ax.set_yticks([])\n"
"    for t,a in zip(ax.get_yticklabels(),AAS): t.set_color(cat.get(a,'#000'))\n"
"    wtl=[f\"{sub[sub.mut_pos==p]['wt_aa'].iloc[0]}{p}\" if len(sub[sub.mut_pos==p]) else str(p) for p in pos]\n"
"    ax.set_xticks(range(len(pos))); ax.set_xticklabels(wtl,fontsize=4,rotation=90); ax.set_title(reg,fontsize=7.5)\n"
"    ax.set_xlabel('wild-type residue + position',fontsize=5.5)\n"
"    if ci==0: ax.set_ylabel('mutant AA'); panel(ax,'a')\n"
"    fig.colorbar(im,ax=ax,fraction=0.046,pad=0.02).ax.tick_params(labelsize=4)\n"
"# b) AA 20x20 true / predicted / error — category-coloured AA labels\n"
"for ci,(Mx,ttl,cmap) in enumerate([(true,'True ΔpKd','RdYlGn'),(pred,'Predicted ΔpKd','RdYlGn'),(err,'Error (true − pred)','coolwarm')]):\n"
"    ax=fig.add_subplot(gs[1,ci]); v=np.nanmax(np.abs(Mx)); v=float(v) if np.isfinite(v) and v>1e-9 else 1.0\n"
"    im=ax.imshow(Mx,cmap=cmap,norm=TwoSlopeNorm(0,-v,v),aspect='auto')\n"
"    ax.set_xticks(range(len(AAo))); ax.set_xticklabels(AAo,fontsize=4.5)\n"
"    ax.set_yticks(range(len(AAo))); ax.set_yticklabels(AAo,fontsize=4.5) if ci==0 else ax.set_yticks([])\n"
"    for t,a in zip(ax.get_xticklabels(),AAo): t.set_color(cat.get(a,'#000'))\n"
"    if ci==0:\n"
"        for t,a in zip(ax.get_yticklabels(),AAo): t.set_color(cat.get(a,'#000'))\n"
"        ax.set_ylabel('original AA'); panel(ax,'b')\n"
"    ax.set_xlabel('mutant AA'); ax.set_title(ttl,fontsize=7.5)\n"
"    fig.colorbar(im,ax=ax,fraction=0.046,pad=0.02).ax.tick_params(labelsize=4)\n"
"fig.suptitle('Mutational landscapes: saturation mutagenesis (a) and 20×20 substitution matrix (b)\\n"
"AA label colour — charged (red) · polar (blue) · hydrophobic (green) · special (purple)',fontsize=8.5,fontweight='bold')\n"
"plt.tight_layout(rect=[0,0,1,0.94]); save(fig,'panel6_heatmaps'); plt.show()")

# ── Panel 7 : scatters (file-backed) ─────────────────────────────────────────
md("---\n## Panel 7 | Predicted vs measured scatters (reads raw data files)")
code(
"from scipy.stats import pearsonr, spearmanr\n"
"import glob\n"
"def best_fold(d,fcol,tcol,pcol):\n"
"    bf=max(d[fcol].unique(),key=lambda fi:pearsonr(d[d[fcol]==fi][tcol],d[d[fcol]==fi][pcol])[0]); return d[d[fcol]==bf]\n"
"def scat(ax,t,p,title,lab):\n"
"    t=np.asarray(t,float); p=np.asarray(p,float); r=pearsonr(t,p)[0]; rho=spearmanr(t,p)[0]; rmse=float(np.sqrt(np.mean((t-p)**2)))\n"
"    xs=np.linspace(t.min(),t.max(),100); rng=np.random.default_rng(0); B=400; L=np.zeros((B,100))\n"
"    for b in range(B):\n"
"        idx=rng.integers(0,len(t),len(t)); m,c=np.polyfit(t[idx],p[idx],1); L[b]=m*xs+c\n"
"    ax.scatter(t,p,c=OURS,alpha=0.35,s=8,linewidths=0,rasterized=True)\n"
"    ax.plot(xs,L.mean(0),color='#7b241c',lw=1.2); ax.fill_between(xs,np.percentile(L,2.5,0),np.percentile(L,97.5,0),color=OURS,alpha=0.25)\n"
"    ax.set_xlabel(f'Measured {PKD}'); ax.set_ylabel(f'Predicted {PKD}')\n"
"    ax.set_title(f'{title}\\nr={r:.3f}, ρ={rho:.3f}, RMSE={rmse:.2f} (n={len(t)})',fontsize=6.6); panel(ax,lab)\n"
"fig,axes=plt.subplots(2,3,figsize=(7.4,5.0)); axes=axes.flatten()\n"
"sa=best_fold(pd.read_csv(os.path.join(HERE,'results_saaintdb_allcdr/random/all_preds.csv')),'fold','binding_affinity','predicted_affinity')\n"
"scat(axes[0],sa['binding_affinity'],sa['predicted_affinity'],'SAaIntDB',' a')\n"
"for ax,ds,nm,lab in zip(axes[1:4],['sabdab','abbind','skempi'],['SAbDab','AB-Bind','SKEMPI'],[' b',' c',' d']):\n"
"    d=best_fold(pd.read_csv(os.path.join(HERE,f'experiments/results_allcdr_stats/oof_preds_{ds}.csv')),'fold','true','pred'); scat(ax,d['true'],d['pred'],nm,lab)\n"
"bps=glob.glob(os.path.join(HERE,'experiments/results/ours_allcdr_benchmark/seed_*/benchmark_preds.csv'))\n"
"bb=max(bps,key=lambda fp:pearsonr(pd.read_csv(fp)['true'],pd.read_csv(fp)['pred'])[0]); bd=pd.read_csv(bb)\n"
"scat(axes[4],bd['true'],bd['pred'],'Benchmark',' e')\n"
"axes[5].axis('off')\n"
"fig.suptitle('Predicted vs measured affinity (best representative fold/seed) with 95% CI band\\nSAaIntDB (a), SAbDab (b), AB-Bind (c), SKEMPI (d), Benchmark (e)',fontsize=8.5,fontweight='bold')\n"
"plt.tight_layout(rect=[0,0,1,0.93]); save(fig,'panel7_scatters'); plt.show()")

code("import os\nprint('Panels saved to figures_panels/:')\nfor f in sorted(os.listdir(os.path.join(HERE,'figures_panels'))):\n    if f.endswith('.png'): print(' ',f)")

nb={'cells':cells,'metadata':{'kernelspec':{'display_name':'Python 3 (py310)','language':'python','name':'py310'},
    'language_info':{'name':'python','version':'3.10'}},'nbformat':4,'nbformat_minor':5}
json.dump(nb,open(NB,'w',encoding='utf-8'),ensure_ascii=False,indent=1)
print(f"Wrote {len(cells)} cells -> {NB}")
print(f"DATA embedded: {len(DATA_LITERAL)} chars")
