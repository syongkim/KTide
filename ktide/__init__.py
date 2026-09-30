"""KTide sequential recursive least squares, Q=0 (paper copy).

Copied from khoa.proj/.../KTide/ktide/. Does not include KHOA operational
helpers (suitability, wet-dry reports). ``fit_ktide_shared`` is the
shared-covariance map-slice update in double precision (v1.1.0).
"""
from .fit import (
    fit_ktide,
    KTideResult,
    fit_ktide_uv,
    KTideVectorResult,
    speed_dir_to_uv,
)
from .constituents import incheon_constituents, ConstituentSet
from .nodal_corrections import compute_nodal
from .benchmark import run_benchmark, build_truth, compare_to_truth
from .grid import fit_ktide_shared

__all__ = [
    "fit_ktide",
    "fit_ktide_uv",
    "KTideResult",
    "KTideVectorResult",
    "speed_dir_to_uv",
    "incheon_constituents",
    "ConstituentSet",
    "compute_nodal",
    "run_benchmark",
    "build_truth",
    "compare_to_truth",
    "fit_ktide_shared",
]
