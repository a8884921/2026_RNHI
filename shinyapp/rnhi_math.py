# -*- coding: utf-8 -*-
"""RNHI math primitives — shared by ALL routes (frozen semantics).
Given: patient composition matrix (n x k closed comps, rows sum to 1) and a center,
score = -Aitchison = -|| CLR(sample) - CLR(center) ||_2  (after multRepl zero handling).
Center is supplied as a composition vector (or CLR vector) computed from a donor table
by comp_center().  Zero handling: multRepl with frac (zCompositions semantics).
"""
import numpy as np
import pandas as pd

LINE12 = ["HSC_MPP", "LMPP", "MyeProg", "Mono", "Gran", "Mega",
          "EarlyEry", "LateEry", "B", "T", "NK", "DC"]
G4 = ["HSC_MPP", "Myeloid", "Erythroid", "Lymphoid"]  # v3 4-lineage order


def closure(x, axis=1):
    x = np.asarray(x, float)
    s = x.sum(axis=axis, keepdims=True)
    if np.any(~np.isfinite(s)) or np.any(s <= 0):
        raise ValueError("closure needs positive row sums")
    return x / s


def multRepl(comp, frac=0.65):
    """zCompositions multRepl(frac) semantics (project-frozen).
    comp: (n,k) closed non-negative composition. Returns replaced comp (rows still sum ~1)."""
    comp = np.asarray(comp, float)
    n, k = comp.shape
    z = comp == 0
    if np.any(comp[~z] <= 0):
        raise ValueError("non-negative required")
    dl = np.full(k, np.nan)
    for j in range(k):
        pos = comp[:, j][comp[:, j] > 0]
        dl[j] = float(np.min(pos)) if pos.size else np.nan
    if np.any(~np.isfinite(dl)):
        raise ValueError("an all-zero column has no detection limit")
    out = comp.copy()
    for i in range(n):
        zi = np.flatnonzero(z[i])
        if zi.size == 0:
            continue
        repl = frac * dl[zi]
        sr = float(repl.sum())
        nzi = np.setdiff1d(np.arange(k), zi)
        sum_nz = float(comp[i, nzi].sum())
        if sum_nz <= 0:
            raise ValueError("row %d fully zero" % i)
        rowsum = float(comp[i].sum())
        scale = (rowsum - sr) / sum_nz
        if scale <= 0 or not np.isfinite(scale):
            raise ValueError("scale<=0 for row %d" % i)
        out[i, zi] = repl
        out[i, nzi] = comp[i, nzi] * scale
    return out


def clr(comp):
    comp = np.asarray(comp, float)
    if np.any(comp <= 0):
        raise ValueError("CLR requires positive components (run multRepl first)")
    lx = np.log(comp)
    return lx - lx.mean(axis=-1, keepdims=True)


def comp_center(donor_wide, method="eq_donor_mean", lineage=LINE12, frac=0.65):
    """donor_wide: DataFrame (donors x lineage cols). Returns CLR-space center vector (k,).
    method:
      'eq_donor_mean' : per-donor composition -> multRepl -> CLR -> arithmetic mean over donors
      'eq_donor_geo'  : geometric mean composition across donors (closure of geometric means) -> CLR
      'pooled'        : pooled counts -> closure -> multRepl -> CLR
      'bary_const'    : donor_wide is a single-row constant composition
    """
    X = donor_wide.reindex(columns=lineage).fillna(0).to_numpy(float)
    # drop donor rows with zero total on the selected axes (no measurable composition in this space)
    X = X[X.sum(axis=1) > 0]
    if len(X) == 0:
        raise ValueError("no donor with positive total on %s" % (lineage,))
    if method == "eq_donor_mean":
        C = closure(X)
        Cr = multRepl(C, frac)
        return clr(Cr).mean(axis=0)
    if method == "eq_donor_geo":
        C = closure(X)
        Cr = multRepl(C, frac)
        gm = np.exp(np.log(Cr).mean(axis=0))
        gm = closure(gm[None, :])[0]
        return clr(gm[None, :])[0]
    if method == "pooled":
        tot = X.sum(axis=0)
        c = closure(tot[None, :])
        cr = multRepl(c, frac)
        return clr(cr)[0]
    if method == "bary_const":
        arr = donor_wide.reindex(columns=lineage).fillna(0).to_numpy(float)
        if arr.ndim == 1:
            arr = arr[None, :]
        c = closure(arr)
        cr = multRepl(c, frac)
        return clr(cr)[0]
    raise ValueError(method)


def rnhi_score(comp, center_clr, frac=0.65, weights=None):
    """comp (n,k closed); center_clr (k,). Return -Aitchison distance (n,).
    weights: optional per-axis weights (1=equal) to down-weight noisy axes."""
    comp = np.asarray(comp, float)
    Cr = multRepl(comp, frac)
    CL = clr(Cr)
    d = CL - center_clr[None, :]
    if weights is not None:
        d = d * np.asarray(weights, float)[None, :]
    return -np.linalg.norm(d, axis=1)


def donor_axis_weights(donor_wide, lineage=LINE12, frac=0.65, power=1.0):
    """Weight = 1/(healthy-donor sd of each CLR axis)^power, standardized so mean=1."""
    X = donor_wide.reindex(columns=lineage).fillna(0).to_numpy(float)
    Cr = multRepl(closure(X), frac)
    CL = clr(Cr)
    sd = CL.std(axis=0)
    sd = np.where(sd < 1e-9, np.nan, sd)
    w = (1.0 / sd) ** power
    w = np.where(np.isnan(w), 1.0, w)
    return w / w.mean()
