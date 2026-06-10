"""
make_architecture_figure.py
===========================
Publication architecture figure for Ours (MutualTriStreamStrong / All-CDR),
drawn to match mutual_strong.py. 4 panels (A-D), 400 DPI, PDF + PNG.
-> figures_publication/fig_architecture.{pdf,png}
"""
import os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'figures_publication'); os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':7,'pdf.fonttype':42,'ps.fonttype':42,
                     'savefig.dpi':400,'figure.dpi':400})

C = {'heavy':'#e74c3c','light':'#f39c12','antigen':'#2ecc71','slate':'#34495e',
     'gate':'#3498db','esm':'#bdc3c7','head':'#8e44ad','bg':'#ffffff','out':'#16a085'}

def box(ax,x,y,w,h,label,fc='#ecf0f1',ec='#2c3e50',fs=6.5,tc='#2c3e50',lw=0.9,style='round',bold=False):
    p=FancyBboxPatch((x,y),w,h,boxstyle=f'{style},pad=0.012,rounding_size=0.02',
                     fc=fc,ec=ec,lw=lw,zorder=2)
    ax.add_patch(p)
    ax.text(x+w/2,y+h/2,label,ha='center',va='center',fontsize=fs,color=tc,zorder=3,
            fontweight='bold' if bold else 'normal',wrap=True)
    return (x+w/2,y+h/2,x,y,w,h)

def arrow(ax,x1,y1,x2,y2,color='#2c3e50',lw=0.9,style='-|>',ls='-'):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle=style,mutation_scale=8,
                 lw=lw,color=color,ls=ls,zorder=1,shrinkA=2,shrinkB=2))

def panel_letter(ax,L):
    ax.text(0.005,0.995,L,transform=ax.transAxes,fontsize=11,fontweight='bold',va='top',ha='left')

fig=plt.figure(figsize=(7.3,8.2))
gs=fig.add_gridspec(2,1,height_ratios=[1.35,1.0],hspace=0.14)
gsb=gs[1].subgridspec(1,3,width_ratios=[1.05,0.95,1.0],wspace=0.13)

# ───────────────────────── Panel A : overall architecture ──────────────────
axA=fig.add_subplot(gs[0]); axA.set_xlim(0,10); axA.set_ylim(0,7.4); axA.axis('off')
panel_letter(axA,'A')
axA.text(5,7.2,'Ours — tri-stream gated cross-attention (final: All-CDR)',
         ha='center',va='top',fontsize=9,fontweight='bold')

# input sequences
def seqbox(name,col,y): return box(axA,0.15,y,1.45,0.85,name,fc=col,ec=col,tc='white',fs=6.6,bold=True)
hY,lY,aY=6.35,4.25,2.0
seqbox('Heavy\nchain',C['heavy'],hY); seqbox('Light\nchain',C['light'],lY); seqbox('Antigen',C['antigen'],aY)

# ESM-2 frozen
for y in (hY,lY,aY):
    box(axA,1.9,y,1.35,0.85,'ESM-2 650M\n(frozen)',fc=C['esm'],ec='#7f8c8d',fs=6.0)
    arrow(axA,1.6,y+0.42,1.9,y+0.42)
# lock icon hint
axA.text(2.57,hY+1.02,'🔒 weights frozen',ha='center',fontsize=5.2,color='#7f8c8d')

# pooling
box(axA,3.5,hY,1.5,0.85,'CDR-aware pool\n(H1+H2+H3)',fc='#fdecea',ec=C['heavy'],fs=5.8,tc=C['heavy'])
box(axA,3.5,lY,1.5,0.85,'mean-pool',fc='#fef5e7',ec=C['light'],fs=6.0,tc='#b9770e')
box(axA,3.5,aY,1.5,0.85,'mean-pool',fc='#eafaf1',ec=C['antigen'],fs=6.0,tc='#1e8449')
for y in (hY,lY,aY): arrow(axA,3.25,y+0.42,3.5,y+0.42)
axA.text(4.25,hY-0.18,'1280-d',ha='center',fontsize=4.8,color='#7f8c8d')

# projections
for y in (hY,lY,aY):
    box(axA,5.25,y,1.45,0.85,'Linear 1280→256\nLN·GELU·Drop',fc='#eaeef2',ec=C['slate'],fs=5.4)
    arrow(axA,5.0,y+0.42,5.25,y+0.42)

# antibody self-attention + fusion gate (heavy+light)
box(axA,7.0,5.3,1.5,0.95,'Heavy↔Light\nself-attention\n(8 heads)',fc='#eef3f8',ec=C['gate'],fs=5.6)
arrow(axA,6.7,hY+0.42,7.0,5.95); arrow(axA,6.7,lY+0.42,7.0,5.55)
box(axA,7.0,4.05,1.5,0.7,'softmax fusion gate',fc='#eef3f8',ec=C['gate'],fs=5.6)
arrow(axA,7.75,5.3,7.75,4.75)
box(axA,7.0,3.15,1.5,0.6,'Antibody repr (256-d)',fc='#f4ecf7',ec=C['head'],fs=5.6,tc=C['head'])
arrow(axA,7.75,4.05,7.75,3.75)
# antigen passes to cross-attn
box(axA,7.0,2.0,1.5,0.6,'Antigen repr (256-d)',fc='#eafaf1',ec=C['antigen'],fs=5.6,tc='#1e8449')
arrow(axA,6.7,aY+0.42,7.0,2.3)

# bidirectional gated cross-attention (x2)
cx=box(axA,8.7,2.55,1.15,1.6,'Bidirectional\nGated\nCross-Attn\n× 2 layers',fc='#d6eaf8',ec=C['gate'],fs=6.0,tc='#21618c',bold=True)
arrow(axA,8.5,3.45,8.7,3.45,color=C['head'])        # antibody -> cross
arrow(axA,8.5,2.3,8.7,2.9,color=C['antigen'])       # antigen -> cross
# loop arrows to indicate bidirectional
axA.text(9.27,4.25,'Ab⇄Ag',ha='center',fontsize=5.2,color='#21618c')

# dual heads + cosine -> pKd
arrow(axA,9.85,3.35,9.85,3.35)  # noop
box(axA,8.7,0.95,1.15,0.8,'Dual heads\nL2-norm',fc='#f4ecf7',ec=C['head'],fs=5.8,tc=C['head'])
arrow(axA,9.27,2.55,9.27,1.75,color=C['gate'])
box(axA,8.7,0.05,1.15,0.65,'cosine → pK$_d$',fc='#d1f2eb',ec=C['out'],fs=6.2,tc='#117864',bold=True)
arrow(axA,9.27,0.95,9.27,0.7,color=C['out'])
axA.text(0.15,0.25,'~2.7M trainable params  (ESM-2 frozen)',fontsize=5.6,color='#7f8c8d',style='italic')

# ───────────────────────── Panel B : gated cross-attention ─────────────────
axB=fig.add_subplot(gsb[0]); axB.set_xlim(0,10); axB.set_ylim(0,10); axB.axis('off')
panel_letter(axB,'B'); axB.text(5,9.6,'Gated Cross-Attention',ha='center',fontsize=8,fontweight='bold')
box(axB,0.6,7.8,3.6,0.9,'Q  (stream 1)',fc='#f4ecf7',ec=C['head'],fs=6.2,tc=C['head'])
box(axB,5.8,7.8,3.6,0.9,'K, V  (stream 2)',fc='#eafaf1',ec=C['antigen'],fs=6.2,tc='#1e8449')
box(axB,1.4,6.0,7.2,0.9,'interaction = tanh(W$_q$Q $\odot$ W$_k$K)',fc='#eef3f8',ec=C['gate'],fs=6.2)
arrow(axB,2.4,7.8,3.6,6.9); arrow(axB,7.6,7.8,6.4,6.9)
box(axB,1.4,4.5,7.2,0.9,'context = interaction $\odot$ W$_v$V',fc='#eef3f8',ec=C['gate'],fs=6.2)
arrow(axB,5.0,6.0,5.0,5.4)
box(axB,0.6,2.9,4.0,0.9,'gate = σ(W$_{gate}$Q)',fc='#d6eaf8',ec=C['gate'],fs=6.2,tc='#21618c',bold=True)
box(axB,5.2,2.9,4.2,0.9,'context $\odot$ gate',fc='#eef3f8',ec=C['gate'],fs=6.2)
arrow(axB,2.4,4.5,2.6,3.8); arrow(axB,6.0,4.5,7.0,3.8); arrow(axB,4.6,3.35,5.2,3.35)
box(axB,1.4,1.1,7.2,0.95,'out = LayerNorm( Q + W$_o$(·) )',fc='#d1f2eb',ec=C['out'],fs=6.4,tc='#117864',bold=True)
arrow(axB,7.3,2.9,5.5,2.05); arrow(axB,2.4,7.8,2.0,2.05,color='#bbb',ls='--')  # residual
axB.text(1.0,2.2,'residual',fontsize=4.8,color='#999',rotation=90)

# ───────────────────────── Panel C : All-CDR pooling ──────────────────────
axC=fig.add_subplot(gsb[1]); axC.set_xlim(0,10); axC.set_ylim(0,10); axC.axis('off')
panel_letter(axC,'C'); axC.text(5,9.6,'CDR-aware pooling',ha='center',fontsize=8,fontweight='bold')
# heavy chain bar with CDR loops
y0=6.8; axC.add_patch(Rectangle((0.6,y0),8.8,0.7,fc='#ecf0f1',ec='#95a5a6',lw=0.8))
for x0,w0,col,lab in [(2.0,1.0,'#f1c40f','H1'),(4.2,1.1,'#2ecc71','H2'),(6.6,1.3,'#9b59b6','H3')]:
    axC.add_patch(Rectangle((x0,y0),w0,0.7,fc=col,ec='none',alpha=0.85))
    axC.text(x0+w0/2,y0+0.35,lab,ha='center',va='center',fontsize=5.6,color='white',fontweight='bold')
axC.text(5.0,y0+1.05,'heavy-chain residues',ha='center',fontsize=5.6,color='#7f8c8d')
# two pooling routes
box(axC,0.6,4.4,4.0,1.0,'mean-pool\nALL residues',fc='#f2f4f4',ec='#95a5a6',fs=5.8)
box(axC,5.4,4.4,4.0,1.0,'mean-pool\nCDR tokens only\n(final model)',fc='#fdecea',ec=C['heavy'],fs=5.6,tc=C['heavy'],bold=True)
arrow(axC,3.0,6.8,2.6,5.4); arrow(axC,7.0,6.8,7.4,5.4,color=C['heavy'])
box(axC,2.6,2.4,4.8,1.0,'heavy embedding (1280-d)',fc='#eaeef2',ec=C['slate'],fs=6.0)
arrow(axC,2.6,4.4,4.4,3.4); arrow(axC,7.4,4.4,5.6,3.4,color=C['heavy'])
axC.text(5.0,1.4,'CDR-aware pooling localizes\nthe paratope signal',ha='center',fontsize=5.6,
         color=C['heavy'],style='italic')

# ───────────────────────── Panel D : objective ────────────────────────────
axD=fig.add_subplot(gsb[2]); axD.set_xlim(0,10); axD.set_ylim(0,10); axD.axis('off')
panel_letter(axD,'D'); axD.text(5,9.6,'Output head & objective',ha='center',fontsize=8,fontweight='bold')
box(axD,1.2,7.6,7.6,0.95,'cosine similarity ∈ [−1, 1]',fc='#d6eaf8',ec=C['gate'],fs=6.2,tc='#21618c')
box(axD,1.2,5.9,7.6,1.0,'pK$_d$ = (cos+1)/2 · (hi−lo) + lo',fc='#d1f2eb',ec=C['out'],fs=6.2,tc='#117864',bold=True)
arrow(axD,5.0,7.6,5.0,6.9)
box(axD,1.2,4.2,7.6,1.0,'MSE( cos , 2(y−lo)/(hi−lo) − 1 )',fc='#fdebd0',ec='#e67e22',fs=6.0,tc='#a04000')
arrow(axD,5.0,5.9,5.0,5.2)
# mini scatter inset
axins=axD.inset_axes([0.18,0.05,0.64,0.30])
import numpy as np
rng=np.random.default_rng(0); t=rng.uniform(5,12,60); p=t+rng.normal(0,0.7,60)
axins.scatter(t,p,s=4,c=C['out'],alpha=0.6,linewidths=0); axins.plot([5,12],[5,12],'k--',lw=0.6)
axins.set_xticks([]); axins.set_yticks([]); axins.set_xlabel('measured',fontsize=4.8); axins.set_ylabel('pred',fontsize=4.8)
for s in axins.spines.values(): s.set_linewidth(0.6)
axD.text(5.0,3.7,'predicted vs measured pK$_d$',ha='center',fontsize=5.0,color='#7f8c8d')

fig.savefig(os.path.join(OUT,'fig_architecture.pdf'),bbox_inches='tight')
fig.savefig(os.path.join(OUT,'fig_architecture.png'),bbox_inches='tight',dpi=400)
print('Saved -> figures_publication/fig_architecture.{pdf,png}')
