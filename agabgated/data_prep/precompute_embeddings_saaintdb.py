#!/usr/bin/env python3
"""
Precompute ESM-2 embeddings for the SAaIntDB dataset.

Reads sequences directly from saaintdb_with_antigen_names.csv (no FASTA needed).
Saves a flat dict  {seq_id: numpy_array}  compatible with CachedEmbeddingLoader.

ID scheme
---------
  heavy_id   : "<PDB_ID>_H_<H_chain_ID>"   e.g. 5zxv_H_H
  light_id   : "<PDB_ID>_L_<L_chain_ID>"   e.g. 5zxv_L_L
  antigen_id : "<PDB_ID>_Ag_<chain(s)>"    e.g. 5zxv_Ag_A

Usage
-----
  python precompute_embeddings_saaintdb.py
  python precompute_embeddings_saaintdb.py --device cpu --batch_size 4
  python precompute_embeddings_saaintdb.py --verify
"""

import os
import pickle
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm

import torch
from transformers import AutoModel, AutoTokenizer


# ---------------------------------------------------------------------------
# ID helpers
# ---------------------------------------------------------------------------

def make_heavy_id(pdb_id, h_chain_id):
    return f"{pdb_id}_H_{h_chain_id}"


def make_light_id(pdb_id, l_chain_id):
    return f"{pdb_id}_L_{l_chain_id}"


def make_antigen_id(pdb_id, ag_chain_ids):
    chains = str(ag_chain_ids).replace(',', '_').replace(' ', '')
    return f"{pdb_id}_Ag_{chains}"


# ---------------------------------------------------------------------------
# ESM-2 embedding computation
# ---------------------------------------------------------------------------

def _compute_batch(seqs, model, tokenizer, device):
    """Mean-pool over non-padding tokens; returns [batch, D] numpy array."""
    with torch.no_grad():
        tokens = tokenizer(
            seqs,
            return_tensors='pt',
            padding=True,
            truncation=True,
            max_length=1024,
        ).to(device)
        out = model(**tokens)
        mask = tokens['attention_mask'].unsqueeze(-1).float()
        emb = (out.last_hidden_state * mask).sum(1) / mask.sum(1)
    return emb.cpu().numpy()


def compute_embeddings_dedup(id_to_seq, model, tokenizer, device, batch_size):
    """
    Compute ESM-2 embeddings for all IDs in id_to_seq.
    Sequences that appear more than once are computed only once.
    Returns dict {id: numpy_array}.
    """
    # Group IDs by identical sequence content (deduplication)
    seq_to_ids: dict[str, list] = {}
    for sid, seq in id_to_seq.items():
        seq_to_ids.setdefault(seq, []).append(sid)

    unique_seqs = list(seq_to_ids.keys())
    n_unique = len(unique_seqs)
    n_total  = len(id_to_seq)
    if n_unique < n_total:
        print(f"    Deduplication: {n_total} IDs → {n_unique} unique sequences")

    # Compute embeddings in batches over unique sequences
    rep_embs: dict[str, np.ndarray] = {}   # seq_content -> embedding
    for start in tqdm(range(0, n_unique, batch_size), desc="    batches"):
        batch_seqs = unique_seqs[start:start + batch_size]
        batch_embs = _compute_batch(batch_seqs, model, tokenizer, device)
        for seq, emb in zip(batch_seqs, batch_embs):
            rep_embs[seq] = emb

    # Map every ID to its embedding (sharing arrays for duplicate sequences)
    result: dict[str, np.ndarray] = {}
    for seq, ids in seq_to_ids.items():
        emb = rep_embs[seq]
        for sid in ids:
            result[sid] = emb

    return result


# ---------------------------------------------------------------------------
# Main precomputation
# ---------------------------------------------------------------------------

def precompute_saaintdb_embeddings(
    csv_path:    str = 'data/saaintdb_with_antigen_names.csv',
    output_path: str = 'data/esm2_embeddings_saaintdb_650M.pkl',
    device:      str = 'cuda',
    batch_size:  int = 8,
    model_name:  str = 'facebook/esm2_t33_650M_UR50D',
):
    print("=" * 70)
    print("PRECOMPUTING ESM-2 EMBEDDINGS — SAaIntDB")
    print("=" * 70)
    print(f"  CSV      : {csv_path}")
    print(f"  Output   : {output_path}")
    print(f"  Model    : {model_name}")
    print(f"  Device   : {device}")
    print(f"  Batch    : {batch_size}")

    df = pd.read_csv(csv_path)
    print(f"\n  Loaded {len(df)} rows")

    # ------------------------------------------------------------------
    # Collect ID → sequence mappings
    # ------------------------------------------------------------------
    heavy_id_to_seq:   dict[str, str] = {}
    light_id_to_seq:   dict[str, str] = {}
    antigen_id_to_seq: dict[str, str] = {}

    skipped = 0
    for _, row in df.iterrows():
        h_id  = make_heavy_id(row['PDB_ID'], row['H_chain_ID'])
        l_id  = make_light_id(row['PDB_ID'], row['L_chain_ID'])
        ag_id = make_antigen_id(row['PDB_ID'], row['Ag_chain_ID(s)'])

        if pd.notna(row['H_seq'])  and str(row['H_seq']).strip():
            heavy_id_to_seq[h_id]   = str(row['H_seq']).strip()
        else:
            skipped += 1

        if pd.notna(row['L_seq'])  and str(row['L_seq']).strip():
            light_id_to_seq[l_id]   = str(row['L_seq']).strip()

        if pd.notna(row['Ag_seq']) and str(row['Ag_seq']).strip():
            antigen_id_to_seq[ag_id] = str(row['Ag_seq']).strip()

    print(f"  Unique heavy    IDs : {len(heavy_id_to_seq)}")
    print(f"  Unique light    IDs : {len(light_id_to_seq)}")
    print(f"  Unique antigen  IDs : {len(antigen_id_to_seq)}")
    if skipped:
        print(f"  Rows with missing H_seq: {skipped}")

    # ------------------------------------------------------------------
    # Load ESM-2
    # ------------------------------------------------------------------
    print(f"\nLoading {model_name} ...")
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    esm_model = AutoModel.from_pretrained(model_name).to(device)
    esm_model.eval()
    print("Model ready.")

    # ------------------------------------------------------------------
    # Compute embeddings
    # ------------------------------------------------------------------
    print("\nHeavy chain embeddings:")
    heavy_embs = compute_embeddings_dedup(
        heavy_id_to_seq, esm_model, tokenizer, device, batch_size)

    print("\nLight chain embeddings:")
    light_embs = compute_embeddings_dedup(
        light_id_to_seq, esm_model, tokenizer, device, batch_size)

    print("\nAntigen embeddings:")
    antigen_embs = compute_embeddings_dedup(
        antigen_id_to_seq, esm_model, tokenizer, device, batch_size)

    # ------------------------------------------------------------------
    # Merge into flat dict and save
    # ------------------------------------------------------------------
    all_embeddings: dict[str, np.ndarray] = {}
    all_embeddings.update(heavy_embs)
    all_embeddings.update(light_embs)
    all_embeddings.update(antigen_embs)

    emb_dim = next(iter(all_embeddings.values())).shape[0]
    print(f"\nTotal embeddings : {len(all_embeddings)}")
    print(f"Embedding dim    : {emb_dim}")

    out_dir = os.path.dirname(output_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(output_path, 'wb') as f:
        pickle.dump(all_embeddings, f)

    size_mb = os.path.getsize(output_path) / (1024 * 1024)
    print(f"Saved to {output_path}  ({size_mb:.1f} MB)")
    return all_embeddings


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------

def verify_cache(
    csv_path:    str = 'data/saaintdb_with_antigen_names.csv',
    cache_path:  str = 'data/esm2_embeddings_saaintdb_650M.pkl',
):
    df = pd.read_csv(csv_path)

    required_heavy   = set()
    required_light   = set()
    required_antigen = set()

    for _, row in df.iterrows():
        required_heavy.add(make_heavy_id(row['PDB_ID'], row['H_chain_ID']))
        required_light.add(make_light_id(row['PDB_ID'], row['L_chain_ID']))
        required_antigen.add(make_antigen_id(row['PDB_ID'], row['Ag_chain_ID(s)']))

    print(f"Loading cache from {cache_path} ...")
    with open(cache_path, 'rb') as f:
        cache = pickle.load(f)

    cached = set(cache.keys())

    miss_h  = required_heavy   - cached
    miss_l  = required_light   - cached
    miss_ag = required_antigen - cached

    print(f"\nRequired heavy   : {len(required_heavy)}  | missing: {len(miss_h)}")
    print(f"Required light   : {len(required_light)}  | missing: {len(miss_l)}")
    print(f"Required antigen : {len(required_antigen)} | missing: {len(miss_ag)}")

    if miss_h:
        print(f"  Sample missing heavy  : {list(miss_h)[:5]}")
    if miss_l:
        print(f"  Sample missing light  : {list(miss_l)[:5]}")
    if miss_ag:
        print(f"  Sample missing antigen: {list(miss_ag)[:5]}")

    ok = not (miss_h or miss_l or miss_ag)
    print("\nCache complete!" if ok else "\nCache INCOMPLETE — re-run without --verify")
    return ok


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Precompute ESM-2 embeddings for SAaIntDB')
    parser.add_argument('--csv',        default='data/saaintdb_with_antigen_names.csv')
    parser.add_argument('--output',     default='data/esm2_embeddings_saaintdb_650M.pkl')
    parser.add_argument('--device',     default='cuda')
    parser.add_argument('--batch_size', type=int, default=8)
    parser.add_argument('--model',      default='facebook/esm2_t33_650M_UR50D')
    parser.add_argument('--verify',     action='store_true',
                        help='Only verify an existing cache, do not recompute')
    args = parser.parse_args()

    if args.verify:
        verify_cache(args.csv, args.output)
    else:
        use_device = args.device if torch.cuda.is_available() else 'cpu'
        precompute_saaintdb_embeddings(
            csv_path    = args.csv,
            output_path = args.output,
            device      = use_device,
            batch_size  = args.batch_size,
            model_name  = args.model,
        )
        verify_cache(args.csv, args.output)
        print("\nDone!  Use the cache in training:")
        print(f"  --cache {args.output}")
