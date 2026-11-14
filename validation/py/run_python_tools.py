#!/usr/bin/env python
"""Run the Python ambient-RNA tools on one exported dataset.

Usage: run_python_tools.py <dataset_dir> [tool ...]     tools: ambidose scar cellbender

<dataset_dir> is written by validation/benchmark_all.R:
  raw/            10x v3 MTX folder: all droplets (cells + empties)
  cells.txt       cell barcodes, in the order the R side expects
  clusters.csv    barcode,cluster
Each tool writes out/<tool>.mtx (genes x cells, the order of raw/features and
cells.txt) and out/<tool>.json (version, seconds, note). A tool that fails
writes only the json, with the error.
"""
import json, sys, time, subprocess, traceback
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp, scipy.io as sio
import scanpy as sc

d = Path(sys.argv[1]); tools = sys.argv[2:] or ["ambidose", "scar", "cellbender"]
out = d / "out"; out.mkdir(exist_ok=True)
cells = [l.strip() for l in open(d / "cells.txt")]
clusters = pd.read_csv(d / "clusters.csv", index_col=0)["cluster"].astype(str)

def load_raw():
    a = sc.read_10x_mtx(d / "raw", var_names="gene_ids", make_unique=False)
    a.X = sp.csr_matrix(a.X)
    return a

def save(tool, X_cells_by_genes, meta):
    sio.mmwrite(str(out / f"{tool}.mtx"), sp.coo_matrix(X_cells_by_genes).T)
    json.dump(meta, open(out / f"{tool}.json", "w"))

def run(tool, fn):
    t0 = time.time()
    try:
        X, meta = fn()
        meta["seconds"] = time.time() - t0
        save(tool, X, meta)
        print(tool, "ok", round(meta["seconds"], 1), "s", flush=True)
    except Exception as e:
        json.dump({"error": f"{type(e).__name__}: {e}", "seconds": time.time() - t0,
                   "trace": traceback.format_exc()[-1500:]}, open(out / f"{tool}.json", "w"))
        print(tool, "FAILED", e, flush=True)

def ambidose_fn():
    import ambidose
    raw = load_raw()
    raw.obs["cluster"] = clusters.reindex(raw.obs_names).values
    filt = raw[cells].copy()
    res = ambidose.denoise(filt, raw=raw, type_key="cluster")
    X = res[cells].X
    rho = None
    for k in ("ambidose_rho", "rho"):
        if k in res.obs: rho = res.obs[k].astype(float).tolist()
    return sp.csr_matrix(X), {"version": ambidose.__version__, "rho": rho,
                              "note": "denoise(filtered, raw=raw, type_key=clusters), defaults"}

def scar_fn():
    import scvi, torch
    from scvi.external import SCAR
    scvi.settings.seed = 0
    raw = load_raw(); filt = raw[cells].copy()
    emp = raw[~raw.obs_names.isin(cells)]
    prof = np.asarray(emp.X.sum(0)).ravel(); prof = prof / prof.sum()
    keep = np.asarray(filt.X.sum(0)).ravel() > 0          # scAR cannot use all-zero genes
    sub = filt[:, keep].copy()
    p = prof[keep]; p = p / p.sum()
    SCAR.setup_anndata(sub)
    model = SCAR(sub, ambient_profile=pd.DataFrame(p, index=sub.var_names, columns=["ambient"]))
    model.train(accelerator="cpu")
    den = sp.csr_matrix(model.get_denoised_counts(n_samples=1).astype(np.float64))
    full = sp.lil_matrix(filt.shape); full[:, np.where(keep)[0]] = den
    return full.tocsr(), {"version": "scvi-tools " + scvi.__version__, "rho": None,
                          "note": "SCAR via scvi-tools, ambient profile from the empties, CPU, default training"}

def cellbender_fn():
    import importlib.metadata as m, h5py
    exe = str(Path(sys.executable).parent / "cellbender")
    h5 = out / "cellbender_out.h5"
    n_total = sum(1 for _ in open(d / "n_droplets.txt")) if False else int(open(d / "n_droplets.txt").read())
    cmd = [exe, "remove-background", "--input", str(d / "raw"), "--output", str(h5),
           "--expected-cells", str(len(cells)), "--total-droplets-included", str(min(n_total, len(cells) + 10000)),
           "--epochs", "150", "--cpu-threads", "4"]
    r = subprocess.run(cmd, cwd=out, capture_output=True, text=True)
    if r.returncode != 0 or not h5.exists():
        raise RuntimeError("cellbender exit %d: %s" % (r.returncode, (r.stderr or r.stdout)[-800:]))
    with h5py.File(h5, "r") as f:
        g = f["matrix"]
        X = sp.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=tuple(g["shape"][:]))
        bc = [b.decode() for b in g["barcodes"][:]]; ids = [b.decode() for b in g["features"]["id"][:]]
    raw_ids = pd.read_csv(d / "raw" / "features.tsv.gz", sep="\t", header=None)[0].tolist()
    X = pd.DataFrame.sparse.from_spmatrix(X, index=ids, columns=bc) if False else X
    gi = pd.Index(ids).get_indexer(raw_ids); ci = pd.Index(bc).get_indexer(cells)
    if (ci < 0).any(): raise RuntimeError("%d cell barcodes missing from CellBender output" % (ci < 0).sum())
    Xc = sp.csr_matrix(X)[gi[gi >= 0]][:, ci]
    full = sp.lil_matrix((len(raw_ids), len(cells))); full[np.where(gi >= 0)[0]] = Xc
    return full.T.tocsr(), {"version": m.version("cellbender"), "rho": None,
                            "note": "remove-background, 150 epochs, CPU, defaults otherwise"}

fns = {"ambidose": ambidose_fn, "scar": scar_fn, "cellbender": cellbender_fn}
for t in tools: run(t, fns[t])
