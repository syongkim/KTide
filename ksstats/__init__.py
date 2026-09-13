"""ksstats sequential mean/var/cov (paper copy).

Copied from /pao1/work/ksstats/src/ksstats/. Formulae: Yoo et al. (2017)
J. Atmos. Oceanic Technol., 34, 2407-2430, Appendix A.
"""
from ksstats.cumstats import (
    CumStats,
    KSMoments,
    KSPairs,
    KSStats,
    KScov,
    KScorr,
    KSmean,
    KSstd,
    KSvar,
    ks_time_mean,
    ks_time_std,
    _CKPT_KEYS,
    batch_numpy_ref,
    batch_same_formula,
    batch_stats,
    compare_batch_vs_cum,
)

__version__ = "0.2.0"

__all__ = [
    "KSMoments",
    "KSPairs",
    "KSStats",
    "CumStats",
    "KSmean",
    "KSstd",
    "KSvar",
    "ks_time_mean",
    "ks_time_std",
    "KScov",
    "KScorr",
    "batch_same_formula",
    "batch_numpy_ref",
    "batch_stats",
    "compare_batch_vs_cum",
    "__version__",
]
