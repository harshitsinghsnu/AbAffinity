"""IG explainability with the FINAL All-CDR SAINTdb model. 400 DPI, separate heatmaps as PDF."""
import os, sys, matplotlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib.pyplot as plt
_orig=plt.savefig
def _sf(*a,**k): k['dpi']=400; return _orig(*a,**k)   # force 400 dpi
plt.savefig=_sf
from AbAffinity.explain import integrated_gradients_explainability_saaintdb as ig3
from AbAffinity.explain import integrated_gradient_twostream as ig2
shared={k:{'light':v['light'],'heavy':v['heavy'],'antigen':v['antigen'],'nanobody':False} for k,v in ig2.COMPLEXES.items()}
ig3.load_complexes=lambda *a,**k: shared
ig3.MODEL_PATH="model_weights/saaintdb_allcdr_random_bestfold.pt"   # final All-CDR (random best fold, val r=0.916)
ig3.OUTPUT_DIR="explainability_results/allcdr_ig"
os.makedirs(ig3.OUTPUT_DIR,exist_ok=True)
ig3.main()
