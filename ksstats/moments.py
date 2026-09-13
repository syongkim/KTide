"""KSMoments — 시간 증분 평균·분산·표준편차 (O(n), 공분산 없음).

문헌: Yoo et al. (2017) JTECH 34, 2407–2430, Appendix A.
KS = Kim–Sequential (시간 증분 + NaN). k-means·K–S test와 무관.

수식 ↔ 이름
    ⟨d⟩     = (1/N_a) Σ d              →  KSmean  /  KSMoments.mean()
    ⟨d²⟩    = [Σ d² − N_a ⟨d⟩²]/(N_a−1) →  KSvar   /  KSMoments.var()
    std(d)  = √⟨d²⟩                    →  KSstd   /  KSMoments.std()

누적량 (체크포인트, 원 시계열 없음)
    c1 = Σ d ,  c2 = Σ d² ,  N = N_a
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import numpy as np

_MOMENT_KEYS = ("c1", "c2", "N")


def n_from_axis(X: np.ndarray, axis: int) -> int:
    X = np.asarray(X)
    if axis == 0:
        return int(X.shape[1])
    if axis == 1:
        return int(X.shape[0])
    raise ValueError("axis must be 0 or 1")


@dataclass
class KSMoments:
    """변수 n개 각각의 시간 증분 평균·분산. 격자 1σ용 (n = ny*nx 가능).

    공분산 배열을 잡지 않음. 페어 통계는 ``KSPairs``.
    """

    n: int
    c1: np.ndarray = field(init=False)
    c2: np.ndarray = field(init=False)
    N: np.ndarray = field(init=False)
    _ckpt_keys: tuple[str, ...] = field(default=_MOMENT_KEYS, init=False, repr=False)

    def __post_init__(self) -> None:
        n = int(self.n)
        if n < 1:
            raise ValueError("n >= 1 required")
        self.n = n
        self.c1 = np.zeros(n, dtype=float)
        self.c2 = np.zeros(n, dtype=float)
        self.N = np.zeros(n, dtype=float)

    def _prepare(self, x: np.ndarray | Iterable[float]) -> tuple[np.ndarray, np.ndarray]:
        d = np.asarray(x, dtype=float).ravel()
        if d.size != self.n:
            raise ValueError(f"expected length {self.n}, got {d.size}")
        finite = np.isfinite(d)
        d0 = np.where(finite, d, 0.0)
        return finite, d0

    def update(self, x: np.ndarray | Iterable[float]) -> None:
        """한 시각 벡터 (길이 n 또는 같은 크기의 2D 장). NaN/Inf는 결측."""
        finite, d0 = self._prepare(x)
        self.N += finite.astype(float)
        self.c1 += d0
        self.c2 += d0 * d0

    def update_many(self, X: np.ndarray, *, axis: int = 1) -> None:
        """여러 시각을 한 번에. axis=1: (n,T) · axis=0: (T,n).

        입력 배열은 호출 측 RAM. 누적기는 여전히 O(n).
        """
        X = np.asarray(X, dtype=float)
        if X.ndim != 2:
            raise ValueError("update_many expects 2D (n,T) or (T,n)")
        if axis == 1:
            if X.shape[0] != self.n:
                raise ValueError(f"expected n={self.n} on axis 0, got {X.shape[0]}")
            finite = np.isfinite(X)
            x0 = np.where(finite, X, 0.0)
            self.N += finite.sum(axis=1).astype(float)
            self.c1 += x0.sum(axis=1)
            self.c2 += (x0 * x0).sum(axis=1)
        elif axis == 0:
            if X.shape[1] != self.n:
                raise ValueError(f"expected n={self.n} on axis 1, got {X.shape[1]}")
            finite = np.isfinite(X)
            x0 = np.where(finite, X, 0.0)
            self.N += finite.sum(axis=0).astype(float)
            self.c1 += x0.sum(axis=0)
            self.c2 += (x0 * x0).sum(axis=0)
        else:
            raise ValueError("axis must be 0 or 1")

    def state_dict(self) -> dict[str, np.ndarray]:
        return {k: getattr(self, k).copy() for k in self._ckpt_keys}

    def load_state_dict(self, state: dict) -> None:
        if "n" in state and int(state["n"]) != self.n:
            raise ValueError(f"n mismatch: state n={state['n']} vs self.n={self.n}")
        for k in self._ckpt_keys:
            if k not in state:
                raise KeyError(f"checkpoint missing {k} (need {type(self).__name__} keys)")
            arr = np.asarray(state[k], dtype=float)
            if arr.shape != getattr(self, k).shape:
                raise ValueError(f"{k} shape {arr.shape} != {getattr(self, k).shape}")
            setattr(self, k, arr.copy())

    def save_checkpoint(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = self.state_dict()
        payload["n"] = np.asarray(self.n, dtype=np.int64)
        payload["kind"] = np.asarray(type(self).__name__)
        np.savez_compressed(path, **payload)
        return path

    @classmethod
    def load_checkpoint(cls, path: str | Path):
        path = Path(path)
        with np.load(path) as z:
            n = int(np.asarray(z["n"]).item())
            obj = cls(n)
            obj.load_state_dict({k: z[k] for k in obj._ckpt_keys} | {"n": n})
        return obj

    def merge(self, other: "KSMoments") -> None:
        if type(other) is not type(self):
            raise TypeError(f"merge same class only, got {type(other).__name__}")
        if other.n != self.n:
            raise ValueError(f"merge n mismatch: {self.n} vs {other.n}")
        for k in self._ckpt_keys:
            setattr(self, k, getattr(self, k) + getattr(other, k))

    def copy(self):
        out = type(self)(self.n)
        out.load_state_dict(self.state_dict() | {"n": self.n})
        return out

    def nbytes(self) -> int:
        return int(sum(getattr(self, k).nbytes for k in self._ckpt_keys))

    def mean(self) -> np.ndarray:
        r"""⟨d⟩ = (1/N_a) Σ d — ``KSmean``."""
        with np.errstate(invalid="ignore", divide="ignore"):
            return self.c1 / self.N

    def var(self) -> np.ndarray:
        r"""⟨d²⟩ = [Σ d² − N_a ⟨d⟩²] / (N_a − 1) — ``KSvar``."""
        mu = self.mean()
        with np.errstate(invalid="ignore", divide="ignore"):
            return self.c2 / (self.N - 1.0) - (self.N / (self.N - 1.0)) * mu**2

    def std(self) -> np.ndarray:
        r"""std(d) = √⟨d²⟩ — ``KSstd``. 격자 1σ."""
        return np.sqrt(self.var())

    def summary(self) -> dict[str, np.ndarray]:
        return {
            "mean": self.mean(),
            "var": self.var(),
            "std": self.std(),
            "N_a": self.N.copy(),
        }


def _fit_moments(X: np.ndarray, axis: int) -> KSMoments:
    X = np.asarray(X, dtype=float)
    ks = KSMoments(n_from_axis(X, axis))
    ks.update_many(X, axis=axis)
    return ks


def KSmean(X: np.ndarray, *, axis: int = 1) -> np.ndarray:
    r"""⟨d⟩ = (1/N_a) Σ d. 내부 ``KSMoments`` (공분산 없음). axis=1 → (n,T)."""
    return _fit_moments(X, axis).mean()


def KSstd(X: np.ndarray, *, axis: int = 1) -> np.ndarray:
    r"""std(d) = √⟨d²⟩, ddof = N_a − 1. 내부 ``KSMoments``."""
    return _fit_moments(X, axis).std()


def KSvar(X: np.ndarray, *, axis: int = 1) -> np.ndarray:
    r"""⟨d²⟩ 표본분산, ddof = N_a − 1. 내부 ``KSMoments``."""
    return _fit_moments(X, axis).var()


def _ks_time_reduce(arr: np.ndarray, *, axis: int, which: str) -> np.ndarray:
    """시간축 축소. ``KSMoments`` only — 격자 (T, ny, nx)에 ``KSPairs`` 금지."""
    a = np.moveaxis(np.asarray(arr, dtype=float), axis, 0)
    spatial = a.shape[1:]
    n = int(np.prod(spatial)) if spatial else 1
    mom = KSMoments(n)
    if spatial:
        mom.update_many(a.reshape(a.shape[0], n), axis=0)
    else:
        mom.update_many(a.reshape(a.shape[0], 1), axis=0)
    out = mom.mean() if which == "mean" else mom.std()
    return out.reshape(spatial) if spatial else out


def ks_time_mean(arr: np.ndarray, *, axis: int = 0) -> np.ndarray:
    """시간 평균. ``np.nanmean(..., axis)`` 과 동일 (유효분 평균). 내부 ``KSMoments``."""
    return _ks_time_reduce(arr, axis=axis, which="mean")


def ks_time_std(arr: np.ndarray, *, axis: int = 0) -> np.ndarray:
    """시간 1σ, ddof = N_a−1. ``np.nanstd(..., axis, ddof=1)`` 과 동일. 격자용."""
    return _ks_time_reduce(arr, axis=axis, which="std")
