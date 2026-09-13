"""ksstats 공개면 — 수식은 README · moments.py · pairs.py 머리글.

    KSMoments  →  KSmean / KSstd / KSvar     (O(n), 격자 1σ)
    KSPairs    →  KScov / KScorr             (O(n²), 정점 페어)
    KSStats    =  KSPairs  (하위 호환 별칭. 격자 전역에 쓰지 말 것)
    CumStats   =  KSPairs
"""
from ksstats.moments import (  # noqa: F401
    KSMoments,
    KSmean,
    KSstd,
    KSvar,
    _MOMENT_KEYS,
    ks_time_mean,
    ks_time_std,
)
from ksstats.pairs import (  # noqa: F401
    KSPairs,
    KScov,
    KScorr,
    _CKPT_KEYS,
    batch_numpy_ref,
    batch_same_formula,
    batch_stats,
    compare_batch_vs_cum,
)

KSStats = KSPairs
CumStats = KSPairs
