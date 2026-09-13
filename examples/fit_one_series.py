"""Minimal KTide call matching Appendix B of the manuscript."""
from __future__ import annotations

from ktide.constituents import incheon_constituents
from ktide.benchmark import build_truth
from ktide.fit import fit_ktide

C = incheon_constituents()
syn = build_truth(C, start="2000-01-01", years=1.0, dt_hours=1.0)
r = fit_ktide(
    syn["time_dnum"],
    syn["zeta_m"],
    C.names,
    C.speed_cpd,
    nodal_mode="amp_phase",
)
i = C.names.index("M2")
print("n_obs", r.n_obs)
print("M2 amp m", float(r.amp_m[i]), "truth", float(C.amp_m[i]))
print("M2 phase deg", float(r.phase_deg[i]), "truth", float(C.phase_deg[i]))
