"""Shared-covariance KTide on a set of cells sampled at the same times (Q = 0).

Cells that share valid sample times, the constituent list, and R carry an
identical parameter covariance P, so one P is kept and every cell keeps its own
coefficient vector m. The identity m = P b (b = sum of G_k^T d_k / R over the
cell's samples) then holds at every cell and the terminal state equals batch
least squares on the same basis.

Precision: P, the time axis, the argument theta, and the basis are float64.
Over a year theta reaches ~1e4 rad, and when the record barely separates two
names (e.g. 365 d with K1-PSI1 in the 53-name list) a single-precision P gives
errors of tens of metres. m may be stored in float32 to save memory.

Missing samples: this function requires every cell to be valid at every time.
Cells with different masks need their own P (see fit_ktide) or grouping by mask;
updating one P whenever some cells are valid breaks m = P b and can diverge.
"""
from __future__ import annotations

import numpy as np

from .nodal_corrections import compute_nodal


def fit_ktide_shared(time_dnum, data, const_names, speed_cpd, *, nodal_mode="amp_phase",
                     R=1e-3, P0=1e5, m_dtype=np.float64):
    """Shared-P recursive least squares over the columns of ``data``.

    Parameters
    ----------
    time_dnum : (nt,) MATLAB datenums (UTC). Phase is local at time_dnum[0].
    data : (nt, n) samples; no NaN allowed.
    const_names, speed_cpd : constituent names and speeds (cycles per day).
    nodal_mode : 'amp_phase' | 'phase_only' | 'none'.
    R, P0 : observation-error variance and initial variance (prior ratio P0/R).
    m_dtype : dtype of the coefficient array (float64 or float32).

    Returns
    -------
    dict with amp (nf, n), phase_deg (nf, n), m (2 nf, n), P (2 nf, 2 nf).
    """
    t = np.asarray(time_dnum, dtype=np.float64).ravel()
    y = np.asarray(data, dtype=np.float64)
    if y.ndim == 1:
        y = y[:, None]
    if not np.isfinite(y).all():
        raise ValueError("fit_ktide_shared needs every cell valid at every time; "
                         "use fit_ktide per cell or group cells by mask")
    names = list(const_names)
    f = np.asarray(speed_cpd, dtype=np.float64).ravel()
    ns = 2 * len(names)
    trel = t - t[0]
    raw = 2.0 * np.pi * np.outer(trel, f)
    mode = nodal_mode.lower()
    if mode == "none":
        gc, gs = np.cos(raw), np.sin(raw)
    else:
        fn, un, _ = compute_nodal(names, t)
        th = raw + un
        if mode == "phase_only":
            gc, gs = np.cos(th), np.sin(th)
        elif mode == "amp_phase":
            gc, gs = fn * np.cos(th), fn * np.sin(th)
        else:
            raise ValueError(f"unknown nodal_mode {nodal_mode!r}")
    P = P0 * np.eye(ns)
    m = np.zeros((ns, y.shape[1]), dtype=m_dtype)
    g = np.empty(ns)
    for k in range(len(t)):
        g[0::2], g[1::2] = gc[k], gs[k]
        pg = P @ g
        K = pg / (g @ pg + R)
        nu = y[k] - g @ m
        m += np.outer(K, nu).astype(m_dtype)
        P -= np.outer(K, pg)
    a, b = m[0::2].astype(np.float64), m[1::2].astype(np.float64)
    return {"amp": np.hypot(a, b), "phase_deg": np.mod(np.degrees(np.arctan2(b, a)), 360.0),
            "m": m, "P": P}
