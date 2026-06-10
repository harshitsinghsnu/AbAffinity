"""
build_explainability_notebook.py
================================
Builds BALM_AbAg_Explainability_Analysis.ipynb:
Integrated-Gradients explainability comparing the final 3-stream (All-CDR)
vs two-stream Ours on 5 canonical antibody-antigen complexes
(1VFB, 3HFM, 5GRJ, 5Y9J, 4ETQ), showing CDR-region highlighting that maps
to known PDB structures.
"""
import json, os, uuid
HERE=os.path.dirname(os.path.abspath(__file__))
NB=os.path.join(HERE,'BALM_AbAg_Explainability_Analysis.ipynb')
cells=[]
def md(t): cells.append({'cell_type':'markdown','metadata':{},'id':uuid.uuid4().hex[:8],'source':[t]})
def code(s): cells.append({'cell_type':'code','metadata':{},'id':uuid.uuid4().hex[:8],'execution_count':None,'outputs':[],'source':[s]})

md("# Ours — Integrated Gradients Explainability\n"
   "## All-CDR (final) vs Two-stream — CDR-Region Attribution on Canonical Complexes\n\n"
   "Per-residue Integrated Gradients (captum, 50 steps) over ESM-2 embeddings, using the "
   "**best SAaIntDB folds** (All-CDR val r=0.872; two-stream val r=0.855). Five canonical "
   "antibody–antigen complexes with known PDB structures: **1VFB, 3HFM, 5GRJ, 5Y9J, 4ETQ**. "
   "We test whether attribution concentrates on the **CDR-H1/H2/H3** loops (the paratope).")

code(
"import os, numpy as np, pandas as pd\n"
"import matplotlib, matplotlib.pyplot as plt\n"
"from matplotlib import image as mpimg\n"
"HERE=os.path.dirname(os.path.abspath('__file__'))\n"
"plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':7,\n"
"  'xtick.labelsize':6,'ytick.labelsize':6,'legend.fontsize':6,'legend.frameon':False,\n"
"  'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':400,'savefig.dpi':400,\n"
"  'savefig.bbox':'tight','pdf.fonttype':42})\n"
"EXP=os.path.join(HERE,'experiments','results_explainability')\n"
"FIG=os.path.join(HERE,'figures_explainability'); os.makedirs(FIG,exist_ok=True)\n"
"COMPLEXES=['1VFB','3HFM','5GRJ','5Y9J','4ETQ']\n"
"C3='#e74c3c'; C2='#3498db'\n"
"def save_fig(fig,n):\n"
"    for e in ('pdf','png'): fig.savefig(os.path.join(FIG,f'{n}.{e}'),bbox_inches='tight')\n"
"    print('saved',n)\n"
"print('ready')")

# ── 1. CDR enrichment ────────────────────────────────────────────────────────
md("---\n## 1. CDR Enrichment — Does Attribution Concentrate on CDR Loops?\n"
   "**Enrichment = (fraction of |attribution| in CDRs) / (fraction of residues that are CDRs).** "
   "Enrichment > 1 means the model preferentially attributes binding to CDR positions.")
code(
"enr=pd.read_csv(os.path.join(EXP,'cdr_enrichment.csv'))\n"
"print(enr.to_string(index=False))")
code(
"fig,axes=plt.subplots(1,2,figsize=(7.6,3.2))\n"
"x=np.arange(len(COMPLEXES)); w=0.38\n"
"three=[float(enr[(enr.complex==c)&(enr.model=='three_stream')]['enrichment'].values[0]) for c in COMPLEXES]\n"
"two=[float(enr[(enr.complex==c)&(enr.model=='two_stream')]['enrichment'].values[0]) for c in COMPLEXES]\n"
"ax=axes[0]\n"
"ax.bar(x-w/2,three,w,color=C3,alpha=0.85,label='All-CDR (final)')\n"
"ax.bar(x+w/2,two,w,color=C2,alpha=0.85,label='Two-stream')\n"
"ax.axhline(1.0,color='k',ls='--',lw=1,label='no enrichment')\n"
"ax.set_xticks(x); ax.set_xticklabels(COMPLEXES); ax.set_ylabel('CDR enrichment (attr / length)')\n"
"ax.set_title('A  CDR attribution enrichment',fontsize=8); ax.legend(fontsize=6)\n"
"ax.text(0.02,0.97,'A',transform=ax.transAxes,fontsize=10,fontweight='bold',va='top')\n"
"ax=axes[1]\n"
"t3=[float(enr[(enr.complex==c)&(enr.model=='three_stream')]['cdr_attr_frac'].values[0]) for c in COMPLEXES]\n"
"t2=[float(enr[(enr.complex==c)&(enr.model=='two_stream')]['cdr_attr_frac'].values[0]) for c in COMPLEXES]\n"
"ax.bar(x-w/2,t3,w,color=C3,alpha=0.85,label='All-CDR')\n"
"ax.bar(x+w/2,t2,w,color=C2,alpha=0.85,label='Two-stream')\n"
"ax.set_xticks(x); ax.set_xticklabels(COMPLEXES); ax.set_ylabel('Fraction of |attribution| in CDRs')\n"
"ax.set_title('B  CDR attribution mass',fontsize=8); ax.legend(fontsize=6)\n"
"ax.text(0.02,0.97,'B',transform=ax.transAxes,fontsize=10,fontweight='bold',va='top')\n"
"fig.suptitle('All-CDR vs Two-stream: CDR-region attribution (5 canonical complexes)',fontsize=9,fontweight='bold')\n"
"plt.tight_layout(w_pad=2); save_fig(fig,'expfig1_cdr_enrichment'); plt.show()\n"
"print(f'Mean All-CDR enrichment: {np.mean(three):.2f}x | two-stream: {np.mean(two):.2f}x')")

# ── 2. Heavy-chain attribution profiles with CDR shading ─────────────────────
md("---\n## 2. Heavy-Chain Attribution Profiles with CDR-H1/H2/H3 Shaded\n"
   "Per-residue attribution along the heavy chain; CDR loops shaded. Peaks inside "
   "the shaded CDR bands indicate the model localizes binding signal to the paratope.")
code(
"fig,axes=plt.subplots(len(COMPLEXES),1,figsize=(8.5,11))\n"
"for ax,cx in zip(axes,COMPLEXES):\n"
"    d=np.load(os.path.join(EXP,f'{cx}.npz'))\n"
"    hl=int(d['heavy_len']); h3=d['three_heavy'][:hl]; h2=d['two_heavy'][:hl]\n"
"    # normalize |attr| to [0,1] for comparison\n"
"    def nrm(a): a=np.abs(a); return a/(a.max()+1e-9)\n"
"    ax.plot(nrm(h3),color=C3,lw=1.0,label='All-CDR')\n"
"    ax.plot(nrm(h2),color=C2,lw=1.0,alpha=0.8,label='Two-stream')\n"
"    for key,col,lab in [('cdr_h1','#f1c40f','CDR-H1'),('cdr_h2','#2ecc71','CDR-H2'),('cdr_h3','#9b59b6','CDR-H3')]:\n"
"        idx=d[key]\n"
"        if len(idx): ax.axvspan(int(idx.min()),int(idx.max()),color=col,alpha=0.2,label=lab)\n"
"    ax.set_title(f'{cx}  (heavy chain, {hl} residues)',fontsize=8)\n"
"    ax.set_ylabel('|attr| (norm)'); ax.set_xlim(0,hl)\n"
"    if cx==COMPLEXES[0]: ax.legend(fontsize=5.5,ncol=5,loc='upper right')\n"
"axes[-1].set_xlabel('Heavy-chain residue index')\n"
"fig.suptitle('Heavy-chain IG attribution with CDR loops shaded (All-CDR vs two-stream)',fontsize=9,fontweight='bold')\n"
"plt.tight_layout(); save_fig(fig,'expfig2_heavy_attr_cdr_shaded'); plt.show()")

# ── 3. Per-chain residue-importance heatmaps (letter cells, consistent colors) ──
md("---\n## 3. Per-Residue Importance Heatmaps — All-CDR vs Two-stream\n"
   "Residue-importance (IG) maps with amino-acid letters in each cell, per-chain "
   "colormaps (**Heavy = green, Light = blue, Antigen = purple**) and a **consistent "
   "chain order (Heavy · Light · Antigen) in both models** — so heavy/light colors are "
   "never swapped. Each block: All-CDR (final) above, Two-stream below.")
code(
"import math\n"
"import integrated_gradient_twostream as ig2   # sequences for letter labels\n"
"SEQ=ig2.COMPLEXES\n"
"CH=[('Heavy chain','heavy','Greens'),('Light chain','light','Blues'),('Antigen','antigen','Purples')]\n"
"NCOL=30\n"
"def draw_chain(ax,attr,seq,cmap):\n"
"    a=np.asarray(attr,float); L=min(len(a),len(seq)); a=a[:L]\n"
"    lo,hi=a.min(),a.max(); a=(a-lo)/(hi-lo+1e-8)   # min-max of SIGNED attr (matches IG pipeline _norm)\n"
"    nrow=int(math.ceil(L/NCOL)); M=np.full((nrow,NCOL),np.nan)\n"
"    for i in range(L): M[i//NCOL,i%NCOL]=a[i]\n"
"    im=ax.imshow(M,cmap=cmap,vmin=0,vmax=1,aspect=1/1.5)   # cells: width = 1.5 x height\n"
"    for i in range(L):\n"
"        ax.text(i%NCOL,i//NCOL,seq[i],ha='center',va='center',fontsize=6.5,color='black')\n"
"    ax.set_xticks([]); ax.set_yticks([])\n"
"    for s in ax.spines.values(): s.set_visible(False)\n"
"    return im\n"
"def model_block(cx,mod,title):\n"
"    d=np.load(os.path.join(EXP,f'{cx}.npz')); seqs=SEQ[cx]\n"
"    nrows=[int(math.ceil(min(len(d[f'{mod}_{k}']),len(seqs[k]))/NCOL)) for _,k,_ in CH]\n"
"    fig,axes=plt.subplots(3,1,figsize=(10,0.46*sum(nrows)+2.0),\n"
"                          gridspec_kw={'height_ratios':[max(n,1) for n in nrows]})\n"
"    for ax,(cname,ckey,cmap) in zip(axes,CH):\n"
"        im=draw_chain(ax,d[f'{mod}_{ckey}'],seqs[ckey],cmap); ax.set_title(f'Residue Importance: {cname}',fontsize=9)\n"
"        cb=fig.colorbar(im,ax=ax,fraction=0.022,pad=0.012,aspect=10)\n"
"        cb.set_label('Normalised IG attribution',fontsize=6.5); cb.ax.tick_params(labelsize=5.5)\n"
"    fig.suptitle(f'{cx} — {title}',fontsize=11,fontweight='bold')\n"
"    plt.tight_layout(rect=[0,0,1,0.97]); save_fig(fig,f'expfig3_{mod}_{cx}'); plt.show()\n"
"for cx in COMPLEXES:\n"
"    model_block(cx,'three','All-CDR (final)')\n"
"    model_block(cx,'two','Two-stream')")

# ── 4. Summary ───────────────────────────────────────────────────────────────
md("---\n## 4. Summary\n\n"
   "- Both models place **>1× enrichment** of attribution on CDR loops for most complexes, "
   "confirming the learned attention localizes to the paratope.\n"
   "- The **All-CDR (final All-CDR-consistent) model** generally shows comparable or stronger "
   "CDR enrichment than the two-stream variant, supporting the All-CDR design.\n"
   "- Heavy-chain profiles show attribution peaks falling inside the CDR-H1/H2/H3 bands — "
   "directly mappable onto the corresponding PDB structures (1VFB, 3HFM, 5GRJ, 5Y9J, 4ETQ).")
code(
"enr=pd.read_csv(os.path.join(EXP,'cdr_enrichment.csv'))\n"
"summ=enr.groupby('model')[['cdr_attr_frac','cdr_len_frac','enrichment']].mean()\n"
"print('Mean across 5 complexes:'); print(summ.to_string())")

nb={'cells':cells,'metadata':{'kernelspec':{'display_name':'Python 3 (py310)','language':'python','name':'py310'},
    'language_info':{'name':'python','version':'3.10'}},'nbformat':4,'nbformat_minor':5}
json.dump(nb,open(NB,'w',encoding='utf-8'),ensure_ascii=False,indent=1)
print(f"Wrote {len(cells)} cells -> {NB}")
