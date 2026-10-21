# -*- coding: utf-8 -*-
"""Python implementation of RNHI (residual normal hematopoiesis integrity).

RNHI of a sample is the negative Aitchison distance between the centered log-ratio
(CLR) transform of its nine-lineage residual normal composition x and the healthy
center c:

    RNHI(x) = -|| clr(x) - c ||_2

Zeros are handled with the multiplicative replacement rule of zCompositions::multRepl:
each zero is replaced with 0.65 times the detection limit of that lineage (its minimum
positive proportion) and the non-zero parts are rescaled so that the composition still
sums to 1. The healthy center is the mean CLR vector of the healthy reference donors,
with zero replacement carried out within the donors.

This module mirrors the R implementation in rnhi_core.R.

Example
-------
>>> import pandas as pd, rnhi_math as rm
>>> tab = pd.read_csv('data/bulk_compendium_composition.csv', index_col='sample_id')
>>> healthy = tab[tab.healthy_reference == 1]
>>> aml = tab[tab.healthy_reference == 0]
>>> center = rm.comp_center(healthy)
>>> dl = rm.detection_limits(rm.closure(aml[rm.LINEAGES].to_numpy(float)))
>>> rnhi = rm.rnhi_score(rm.closure(aml[rm.LINEAGES].to_numpy(float)), center, dl=dl)
"""
import numpy as np
import pandas as pd

LINEAGES = ["MyeProg", "Mono", "Gran", "Mega", "EarlyEry", "LateEry", "B", "NK", "DC"]
LINE9 = LINEAGES
ZERO_REPLACEMENT_FACTOR = 0.65


def closure(x):
    """Rescale each row of a non-negative matrix to sum to 1."""
    x = np.asarray(x, float)
    if x.ndim == 1:
        x = x[None, :]
    if np.any(~np.isfinite(x)) or np.any(x < 0):
        raise ValueError("abundances must be finite and non-negative")
    s = x.sum(axis=1, keepdims=True)
    if np.any(s <= 0):
        raise ValueError("at least one sample has no detected lineage; RNHI is undefined")
    return x / s


def detection_limits(comp):
    """Minimum positive proportion of each lineage (NaN if the lineage is never detected)."""
    comp = np.asarray(comp, float)
    dl = np.full(comp.shape[1], np.nan)
    for j in range(comp.shape[1]):
        pos = comp[:, j][comp[:, j] > 0]
        if pos.size:
            dl[j] = float(pos.min())
    return dl


def multRepl(comp, frac=ZERO_REPLACEMENT_FACTOR, dl=None):
    """Multiplicative zero replacement of closed compositions (rows sum to 1).

    If dl is None, detection limits are estimated from comp itself; otherwise the
    supplied limits (e.g. the fixed limits of the bulk AML compendium) are used.
    """
    if not 0 < frac < 1:
        raise ValueError("the zero-replacement factor must lie between 0 and 1")
    comp = np.asarray(comp, float)
    dl = detection_limits(comp) if dl is None else np.asarray(dl, float)
    zero = comp == 0
    needs = zero.any(axis=0)
    if np.any(needs & ~(np.isfinite(dl) & (dl > 0))):
        raise ValueError("no detection limit for a lineage that is zero in every sample of the group")
    out = comp.copy()
    for i in range(comp.shape[0]):
        zi = np.flatnonzero(zero[i])
        if zi.size == 0:
            continue
        nzi = np.flatnonzero(~zero[i])
        repl = frac * dl[zi]
        scale = (comp[i].sum() - repl.sum()) / comp[i, nzi].sum()
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError("zero replacement leaves no positive mass for row %d" % i)
        out[i, zi] = repl
        out[i, nzi] = comp[i, nzi] * scale
    return out


def clr(comp):
    """Centered log-ratio transform of strictly positive compositions."""
    comp = np.asarray(comp, float)
    if np.any(comp <= 0):
        raise ValueError("CLR requires positive parts (apply multRepl first)")
    lx = np.log(comp)
    return lx - lx.mean(axis=-1, keepdims=True)


def comp_center(donor_wide, method="eq_donor_mean", lineage=LINEAGES, frac=ZERO_REPLACEMENT_FACTOR):
    """Healthy center: mean CLR vector of the healthy donors (donors weighted equally).

    donor_wide: DataFrame with one row per donor and the lineage columns.
    Zero replacement is carried out within the donors.
    """
    if method != "eq_donor_mean":
        raise ValueError("the RNHI healthy center is the equal-donor mean CLR vector ('eq_donor_mean')")
    missing = [c for c in lineage if c not in donor_wide.columns]
    if missing:
        raise ValueError("missing lineage columns: %s" % ", ".join(missing))
    X = donor_wide[list(lineage)].to_numpy(float)
    return clr(multRepl(closure(X), frac)).mean(axis=0)


healthy_center = comp_center


def rnhi_score(comp, center_clr, frac=ZERO_REPLACEMENT_FACTOR, dl=None):
    """RNHI (negative Aitchison distance to the healthy center) of closed compositions."""
    d = clr(multRepl(comp, frac, dl)) - np.asarray(center_clr, float)[None, :]
    return -np.sqrt((d ** 2).sum(axis=1))


def score_table(table, healthy=None, dl=None, frac=ZERO_REPLACEMENT_FACTOR):
    """Score a wide composition table (DataFrame indexed by sample_id with the nine lineage columns).

    healthy: DataFrame of healthy reference donors (same columns) used for the center.
    dl: detection limits for the query samples; None estimates them from `table`.
    Returns a Series of RNHI values indexed like `table`.
    """
    if healthy is None:
        raise ValueError("a healthy reference table is required")
    center = comp_center(healthy, lineage=LINEAGES, frac=frac)
    comp = closure(table[LINEAGES].to_numpy(float))
    return pd.Series(rnhi_score(comp, center, frac, dl), index=table.index, name="RNHI")
