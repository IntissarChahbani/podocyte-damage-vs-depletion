import numpy as np
import pandas as pd


def auc_score(expr, genes, auc_max_frac=0.05):
    """AUC of a gene set per sample. expr: genes x samples (index = gene IDs)."""
    ranks = expr.rank(axis=0, ascending=False, method="average")   # 1 = highest
    max_rank = int(np.ceil(auc_max_frac * expr.shape[0]))
    sub = ranks.loc[genes]
    contrib = (max_rank - sub).clip(lower=0)                        # only genes inside the top max_rank count
    return contrib.sum(axis=0) / (len(genes) * max_rank)


def pds_score(expr, up, down, scale=True):
    """Podocyte damage score: HIGHER = MORE damage.
    DS = (n_up/N)*AUC_up - (n_down/N)*AUC_down, N = n_up + n_down (genes present)."""
    up = [g for g in up if g in expr.index]
    down = [g for g in down if g in expr.index]
    n_up, n_down = len(up), len(down)
    n_total = n_up + n_down
    ds_up = auc_score(expr, up) if n_up else 0.0
    ds_down = auc_score(expr, down) if n_down else 0.0
    ds = (n_up / n_total) * ds_up - (n_down / n_total) * ds_down
    if scale:
        ds = (ds - ds.mean()) / ds.std(ddof=1)
    return ds
