"""Builds BALM_AbAg_Transfer_AllCDR.ipynb: zero/few-shot transfer of the final
All-CDR SAaIntDB model onto 6 external datasets (3 seeds, 10/20/30%). 400 DPI PDF."""
import json, os, uuid
HERE=os.path.dirname(os.path.abspath(__file__))
NB=os.path.join(HERE,'BALM_AbAg_Transfer_AllCDR.ipynb')
cells=[]
def md(t): cells.append({'cell_type':'markdown','metadata':{},'id':uuid.uuid4().hex[:8],'source':[t]})
def code(s): cells.append({'cell_type':'code','metadata':{},'id':uuid.uuid4().hex[:8],'execution_count':None,'outputs':[],'source':[s]})
md("# Ours (final All-CDR) — Zero-shot & Few-shot Transfer\n"
   "## Six external single-antigen datasets: 1mlc, 1n8z, 4fqi, koenig, trastuzumab, warszawski\n\n"
   "All-CDR SAaIntDB weights. Zero-shot = 10-fold All-CDR ensemble; few-shot = fine-tune at "
   "**10/20/30%** with **3 seeds** (mean ± std). 4fqi subsampled to 4000 variants.")
code(
"import os, numpy as np, pandas as pd, matplotlib.pyplot as plt\n"
"HERE=os.path.dirname(os.path.abspath('__file__'))\n"
"plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':7,\n"
"  'xtick.labelsize':6,'ytick.labelsize':6,'legend.fontsize':6,'legend.frameon':False,\n"
"  'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':400,'savefig.dpi':400,\n"
"  'savefig.bbox':'tight','pdf.fonttype':42})\n"
"FIG=os.path.join(HERE,'figures_transfer'); os.makedirs(FIG,exist_ok=True)\n"
"def save_fig(f,n):\n"
"    for e in ('pdf','png'): f.savefig(os.path.join(FIG,f'{n}.{e}'),bbox_inches='tight')\n"
"    print('saved',n)\n"
"agg=pd.read_csv(os.path.join(HERE,'experiments','results_transfer_allcdr','aggregated.csv'))\n"
"DS=['1mlc','1n8z','4fqi','koenig','trastuzumab','warszawski']\n"
"print(agg.to_string(index=False))")
md("---\n## 1. Few-shot Learning Curves (per dataset, 3 seeds)")
code(
"fig,axes=plt.subplots(2,3,figsize=(10,6)); axes=axes.flatten()\n"
"for ax,ds,letter in zip(axes,DS,'ABCDEF'):\n"
"    a=agg[agg.dataset==ds].sort_values('fraction')\n"
"    fs=a[a.fraction>0]\n"
"    zs=a[a.fraction==0]['pearson_mean'].values[0]\n"
"    ax.axhline(zs,color='#95a5a6',ls=':',lw=1,label=f'zero-shot ({zs:.2f})')\n"
"    ax.errorbar(fs.fraction*100,fs.pearson_mean,yerr=fs.pearson_std,fmt='o-',color='#e74c3c',lw=1.8,ms=5,capsize=3,label='few-shot')\n"
"    ax.set_xlabel('Training data (%)'); ax.set_ylabel('Pearson r' if letter in 'AD' else '')\n"
"    ax.set_title(f'{letter}  {ds}',fontsize=8); ax.legend(fontsize=6); ax.set_xticks([10,20,30])\n"
"    ax.text(0.02,0.97,letter,transform=ax.transAxes,fontsize=10,fontweight='bold',va='top')\n"
"fig.suptitle('All-CDR transfer: few-shot learning curves (mean ± std, 3 seeds)',fontsize=9,fontweight='bold')\n"
"plt.tight_layout(); save_fig(fig,'transfig1_fewshot_curves'); plt.show()")
md("---\n## 2. Zero-shot vs Few-shot Summary (all metrics)")
code(
"fig,axes=plt.subplots(1,3,figsize=(10,3.2))\n"
"x=np.arange(len(DS)); w=0.2\n"
"fracs=[0.0,0.1,0.2,0.3]; cols=['#95a5a6','#f39c12','#3498db','#e74c3c']; labs=['zero','10%','20%','30%']\n"
"for ax,metric,ylab,letter in [(axes[0],'pearson_mean','Pearson r','A'),(axes[1],'spearman_mean','Spearman \u03c1','B'),(axes[2],'rmse_mean','RMSE','C')]:\n"
"    for j,(fr,c,lb) in enumerate(zip(fracs,cols,labs)):\n"
"        vals=[float(agg[(agg.dataset==d)&(agg.fraction==fr)][metric].values[0]) if len(agg[(agg.dataset==d)&(agg.fraction==fr)]) else np.nan for d in DS]\n"
"        ax.bar(x+(j-1.5)*w,vals,w,color=c,alpha=0.85,label=lb if letter=='A' else None)\n"
"    ax.set_xticks(x); ax.set_xticklabels(DS,rotation=30,ha='right',fontsize=6); ax.set_ylabel(ylab)\n"
"    ax.set_title(f'{letter}  {ylab}',fontsize=8)\n"
"    if metric=='rmse_mean': ax.set_yscale('log')\n"
"    if letter=='A': ax.legend(fontsize=6,title='train frac')\n"
"    ax.text(0.02,0.97,letter,transform=ax.transAxes,fontsize=10,fontweight='bold',va='top')\n"
"fig.suptitle('All-CDR transfer across 6 datasets — zero-shot vs few-shot',fontsize=9,fontweight='bold')\n"
"plt.tight_layout(); save_fig(fig,'transfig2_summary_allmetrics'); plt.show()")
md("---\n## 3. Summary Table")
code("print(agg.round(4).to_string(index=False))")
nb={'cells':cells,'metadata':{'kernelspec':{'display_name':'Python 3 (py310)','language':'python','name':'py310'},'language_info':{'name':'python','version':'3.10'}},'nbformat':4,'nbformat_minor':5}
json.dump(nb,open(NB,'w',encoding='utf-8'),ensure_ascii=False,indent=1)
print(f"Wrote {len(cells)} cells -> {NB}")
