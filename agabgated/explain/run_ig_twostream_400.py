"""Two-stream IG at 400 DPI, separate heatmaps as PDF (ConcatTwoStream, best SAINTdb fold)."""
import os, sys, torch, matplotlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib.pyplot as plt
_orig=plt.savefig
def _sf(*a,**k): k['dpi']=400; return _orig(*a,**k)
plt.savefig=_sf
from agabgated.explain import integrated_gradient_twostream as ig
from agabgated.models.two_stream_mutualstrong import ConcatTwoStream
def load_concat(p,device=ig.DEVICE):
    ck=torch.load(p,map_location=device); cfg=ck.get('config',{})
    m=ConcatTwoStream(esm_dim=1280,projected_size=cfg.get('projected_size',256),num_heads=cfg.get('num_heads',8),
        dropout=0.0,n_layers=cfg.get('n_layers',2),device=device).to(device)
    m.load_state_dict(ck['model_state_dict']); m.eval()
    return m, ck.get('pkd_bounds',(5.0,12.0))
ig.load_model=load_concat
ig.MODEL_PATH="results_twostream_saaintdb/fold_01/model.pt"
ig.OUTPUT_DIR="explainability_results/two_stream_ig"
os.makedirs(ig.OUTPUT_DIR,exist_ok=True)
ig.main()
