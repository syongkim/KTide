"""Resume KTide from the parameter checkpoint (m, P) plus calendar epoch.

Matches section 2 of the manuscript: twelve 30-day increments equal one
360-day pass when Q = 0. The pair (m, P) is the sufficient statistic;
m alone cannot be resumed. ksstats save_checkpoint / load_checkpoint store
moment sums, not this harmonic state.
"""
from __future__ import annotations

import numpy as np

from ktide.constituents import incheon_constituents
from ktide.benchmark import build_truth
from ktide.fit import fit_ktide

C = incheon_constituents()
syn = build_truth(C, start="2000-01-01", years=1.0, dt_hours=1.0)
n = 12 * 30 * 24  # twelve 30-day hourly increments
dnum = np.asarray(syn["time_dnum"], float)[:n]
y = np.asarray(syn["zeta_m"], float)[:n]
epoch = float(dnum[0])

one = fit_ktide(
    dnum,
    y,
    C.names,
    C.speed_cpd,
    nodal_mode="amp_phase",
    epoch_dnum=epoch,
)

n_chunk = 30 * 24
m = None
P = None
r = None
for i in range(12):
    sl = slice(i * n_chunk, (i + 1) * n_chunk)
    r = fit_ktide(
        dnum[sl],
        y[sl],
        C.names,
        C.speed_cpd,
        nodal_mode="amp_phase",
        epoch_dnum=epoch,
        m_init=m,
        P_init=P,
    )
    m = r.model_vector
    P = r.posterior_covariance

dm = float(np.max(np.abs(one.model_vector - r.model_vector)))
dP = float(np.max(np.abs(one.posterior_covariance - r.posterior_covariance)))
print("max |Δm|", dm)
print("max |ΔP|", dP)
print("identity", dm < 1e-8 and dP < 1e-8)
