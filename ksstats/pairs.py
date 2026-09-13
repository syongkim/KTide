"""KSPairs — 시간 증분 공분산·상관 (O(n²) 페어 누적).

문헌: Yoo et al. (2017) JTECH 34, 2407–2430, Appendix A.

수식 ↔ 이름
    ⟨d_p d_q⟩  (분모 N_c−1, 주변 평균 ⟨d_p⟩,⟨d_q⟩)  →  KScov  /  KSPairs.cov()
    ρ_pq = ⟨d_p d_q⟩ / (σ_p σ_q)                   →  KScorr /  KSPairs.corr()

평균·std는 부모 ``KSMoments`` (주변 평균이 공분산에 필요).

격자 전역(n = ny*nx)은 \(O(n^2)\). RAM·목적 확인 후 사용자 판단. 1σ만이면 `KSMoments`.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

import numpy as np

from ksstats.moments import KSMoments, n_from_axis

_CKPT_KEYS = ("c1", "c2", "N", "c3", "c4", "c5", "Nc")


def moments_checkpoint_bytes(n: int) -> int:
    """``KSMoments`` 체크포인트 이론 RAM (float64, c1·c2·N)."""
    return 3 * int(n) * 8


def pairs_checkpoint_bytes(n: int) -> int:
    """``KSPairs`` 체크포인트 이론 RAM (float64). ``nbytes()`` 와 동일."""
    n = int(n)
    npairs = n * (n + 1) // 2
    return 3 * n * 8 + 4 * npairs * 8


def _pairs_theory_bytes(n: int) -> int:
    return pairs_checkpoint_bytes(n)


def _combsrep(n: int) -> tuple[np.ndarray, np.ndarray]:
    """반복 허용 2-조합 (하삼각 포함 대각). n=2 → (0,0),(0,1),(1,1)."""
    g1, g2 = [], []
    for i in range(n):
        for j in range(i, n):
            g1.append(i)
            g2.append(j)
    return np.asarray(g1, dtype=int), np.asarray(g2, dtype=int)


def _vech_to_matrix(vech: np.ndarray, n: int) -> np.ndarray:
    out = np.full((n, n), np.nan, dtype=float)
    k = 0
    for i in range(n):
        for j in range(i, n):
            out[i, j] = out[j, i] = vech[k]
            k += 1
    return out


@dataclass
class KSPairs(KSMoments):
    """n개 변수 사이 모든 페어의 시간 증분 공분산.

    누적량: KSMoments의 c1,c2,N  +  c3=Σ d_p d_q, c4,c5, Nc=N_c.
    """

    c3: np.ndarray = field(init=False)
    c4: np.ndarray = field(init=False)
    c5: np.ndarray = field(init=False)
    Nc: np.ndarray = field(init=False)
    g1: np.ndarray = field(init=False, repr=False)
    g2: np.ndarray = field(init=False, repr=False)

    def __post_init__(self) -> None:
        super().__post_init__()
        if self.n >= 512:
            import warnings

            warnings.warn(
                f"KSPairs(n={self.n}) checkpoint ≈ {pairs_checkpoint_bytes(self.n) / 1024**3:.2f} GiB "
                f"(O(n²); see ksstats README · memory). Grid 1σ only needs KSMoments (~{moments_checkpoint_bytes(self.n) / 1024**2:.2f} MiB).",
                UserWarning,
                stacklevel=2,
            )
        self._ckpt_keys = _CKPT_KEYS
        npairs = self.n * (self.n + 1) // 2
        self.c3 = np.zeros(npairs, dtype=float)
        self.c4 = np.zeros(npairs, dtype=float)
        self.c5 = np.zeros(npairs, dtype=float)
        self.Nc = np.zeros(npairs, dtype=float)
        self.g1, self.g2 = _combsrep(self.n)

    def update(self, x: np.ndarray | Iterable[float]) -> None:
        finite, d0 = self._prepare(x)
        self.N += finite.astype(float)
        self.c1 += d0
        self.c2 += d0 * d0
        both = finite[self.g1] & finite[self.g2]
        self.Nc += both.astype(float)
        self.c3 += d0[self.g1] * d0[self.g2]
        kk = np.where(both)[0]
        if kk.size:
            self.c4[kk] += d0[self.g1[kk]]
            self.c5[kk] += d0[self.g2[kk]]

    def update_many(self, X: np.ndarray, *, axis: int = 1) -> None:
        """페어 누적은 시각 루프. (부모의 벡터 적률 경로는 c3–Nc를 안 갱신함.)"""
        X = np.asarray(X, dtype=float)
        if axis == 0:
            for t in range(X.shape[0]):
                self.update(X[t])
        elif axis == 1:
            for t in range(X.shape[1]):
                self.update(X[:, t])
        else:
            raise ValueError("axis must be 0 or 1")

    def cov_vech(self) -> np.ndarray:
        r"""⟨d_p d_q⟩ 표본공분산, 분모 (N_c − 1), 평균은 주변 ⟨d_p⟩, ⟨d_q⟩."""
        mu = self.mean()
        g1, g2 = self.g1, self.g2
        with np.errstate(invalid="ignore", divide="ignore"):
            return (
                self.c3 / (self.Nc - 1.0)
                - mu[g1] * self.c5 / (self.Nc - 1.0)
                - mu[g2] * self.c4 / (self.Nc - 1.0)
                + mu[g1] * mu[g2] * self.Nc / (self.Nc - 1.0)
            )

    def cov(self) -> np.ndarray:
        r"""공분산 행렬 — ``KScov``."""
        return _vech_to_matrix(self.cov_vech(), self.n)

    def corr_vech(self) -> np.ndarray:
        v = self.var()
        c = self.cov_vech()
        with np.errstate(invalid="ignore", divide="ignore"):
            return c / np.sqrt(v[self.g1] * v[self.g2])

    def corr(self) -> np.ndarray:
        r"""ρ_pq = ⟨d_p d_q⟩ / (σ_p σ_q) — ``KScorr``."""
        return _vech_to_matrix(self.corr_vech(), self.n)

    def summary(self) -> dict[str, np.ndarray]:
        out = super().summary()
        out["cov"] = self.cov()
        out["corr"] = self.corr()
        out["N_c_vech"] = self.Nc.copy()
        return out


def _fit_pairs(X: np.ndarray, axis: int) -> KSPairs:
    X = np.asarray(X, dtype=float)
    ks = KSPairs(n_from_axis(X, axis))
    ks.update_many(X, axis=axis)
    return ks


def KScov(X: np.ndarray, *, axis: int = 1) -> np.ndarray:
    r"""⟨d_p d_q⟩. 내부 ``KSPairs``. 소수 변수만."""
    return _fit_pairs(X, axis).cov()


def KScorr(X: np.ndarray, *, axis: int = 1) -> np.ndarray:
    r"""ρ_pq = KScov / (KSstd_p KSstd_q). 내부 ``KSPairs``."""
    return _fit_pairs(X, axis).corr()


def batch_same_formula(X: np.ndarray, *, axis: int = 1) -> dict[str, np.ndarray]:
    """전체 시계열을 메모리에 올린 뒤 동일 충분통계로 한 번에 계산."""
    X = np.asarray(X, dtype=float)
    if axis == 0:
        X = X.T
    elif axis != 1:
        raise ValueError("axis must be 0 or 1")
    n = X.shape[0]
    cum = KSPairs(n)
    finite = np.isfinite(X)
    cum.N = finite.sum(axis=1).astype(float)
    X0 = np.where(finite, X, 0.0)
    cum.c1 = X0.sum(axis=1)
    cum.c2 = (X0 * X0).sum(axis=1)
    both = finite[cum.g1] & finite[cum.g2]
    cum.Nc = both.sum(axis=1).astype(float)
    cum.c3 = (X0[cum.g1] * X0[cum.g2]).sum(axis=1)
    mask = both.astype(float)
    cum.c4 = (X0[cum.g1] * mask).sum(axis=1)
    cum.c5 = (X0[cum.g2] * mask).sum(axis=1)
    return cum.summary()


def batch_numpy_ref(X: np.ndarray, *, axis: int = 1) -> dict[str, np.ndarray]:
    """참고: nanmean / nanvar(ddof=1) / 페어 평균 기반 nancov."""
    X = np.asarray(X, dtype=float)
    if axis == 1:
        A = X.T
    elif axis == 0:
        A = X
    else:
        raise ValueError("axis must be 0 or 1")

    n = A.shape[1]
    mean = np.nanmean(A, axis=0)
    var = np.nanvar(A, axis=0, ddof=1)
    cov = np.full((n, n), np.nan, dtype=float)
    for i in range(n):
        for j in range(i, n):
            m = np.isfinite(A[:, i]) & np.isfinite(A[:, j])
            ni = int(m.sum())
            if ni < 2:
                continue
            xi, xj = A[m, i], A[m, j]
            cov[i, j] = cov[j, i] = np.sum((xi - xi.mean()) * (xj - xj.mean())) / (ni - 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        s = np.sqrt(np.diag(cov))
        corr = cov / np.outer(s, s)
    return {"mean": mean, "var": var, "std": np.sqrt(var), "cov": cov, "corr": corr}


batch_stats = batch_same_formula


def compare_batch_vs_cum(
    X: np.ndarray, *, axis: int = 1, rtol: float = 1e-10, atol: float = 1e-10
) -> dict:
    batch = batch_same_formula(X, axis=axis)
    n = batch["mean"].size
    cum = KSPairs(n)
    cum.update_many(X, axis=axis)
    online = cum.summary()
    ref = batch_numpy_ref(X, axis=axis)

    diffs = {}
    ok = True
    for key in ("mean", "var", "cov", "corr"):
        d = np.nanmax(np.abs(batch[key] - online[key]))
        diffs[key] = float(d) if np.isfinite(d) else float("nan")
        if not np.allclose(batch[key], online[key], rtol=rtol, atol=atol, equal_nan=True):
            ok = False

    diffs_vs_numpy = {}
    for key in ("mean", "var", "cov", "corr"):
        d = np.nanmax(np.abs(online[key] - ref[key]))
        diffs_vs_numpy[key] = float(d) if np.isfinite(d) else float("nan")

    return {
        "ok": ok,
        "max_abs_diff": diffs,
        "max_abs_diff_vs_numpy_nancov": diffs_vs_numpy,
        "batch": batch,
        "cum": online,
        "numpy_ref": ref,
    }
