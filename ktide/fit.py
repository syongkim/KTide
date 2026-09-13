"""
fit.py
Sequential (recursive least-squares / Kalman-filter) harmonic estimator.
Ported from khoa_fit_ktide.m
"""
from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional

from .nodal_corrections import compute_nodal


@dataclass
class KTideResult:
    names: List[str]
    speed_cpd: np.ndarray
    nodal_mode: str
    amp_m: np.ndarray
    phase_deg: np.ndarray
    model_vector: np.ndarray        # [cos1, sin1, cos2, sin2, …]
    posterior_covariance: np.ndarray
    time_hist: np.ndarray
    amp_hist: np.ndarray            # (n_store, nf)
    phase_hist_deg: np.ndarray      # (n_store, nf)
    node_longitude_deg: np.ndarray  # (nt,) lunar node longitude
    n_obs: int = 0                  # 실제 사용된 관측 수 (NaN 제외) = data_length
    gap_fraction: float = 0.0       # 결측 비율 (0=완전, 1=전부 결측)
    data_length: int = 0            # 내부변수: 자료 있으면 +1 누적 (n_obs와 동일)
    n_kalman_steps: int = 0         # 순차 필터가 지난 시간 스텝 수 (= nt)
    data_fraction: float = 1.0      # data_length / n_kalman_steps


def _cart_to_amp_phase(m: np.ndarray):
    """Convert interleaved [cos, sin, cos, sin, …] vector to amp/phase."""
    nf = len(m) // 2
    amp = np.empty(nf)
    phase = np.empty(nf)
    for j in range(nf):
        c = m[2 * j]
        s = m[2 * j + 1]
        amp[j] = np.hypot(c, s)
        phase[j] = np.mod(np.degrees(np.arctan2(s, c)), 360.0)
    return amp, phase


def fit_ktide(
    time_dnum: np.ndarray,
    zeta_m: np.ndarray,
    const_names: List[str],
    speed_cpd: np.ndarray,
    *,
    nodal_mode: str = 'amp_phase',
    R: float = 1e-3,
    P0: float = 1e5,
    store_stride: Optional[int] = None,
    epoch_dnum: Optional[float] = None,
    m_init: Optional[np.ndarray] = None,
    P_init: Optional[np.ndarray] = None,
) -> KTideResult:
    """
    Generalised KTide harmonic estimator (sequential Kalman update).

    NaN 값은 보간 없이 자동으로 건너뜁니다.
    결측이 있어도 별도 전처리 없이 그대로 입력하면 됩니다.
    결측 스텝에서는 P도 업데이트되지 않아 불확실도가 정직하게 유지됩니다.

    Parameters
    ----------
    time_dnum : ndarray, shape (nt,)
        MATLAB serial datenums. 등간격일 필요 없음.
        관측 시각만 정확히 제공하면 됩니다.
    zeta_m : ndarray, shape (nt,)
        해면 수위 또는 임의의 스칼라 관측 [m].
        결측은 NaN으로 표시. 보간 불필요.
    const_names : list of str
        분조 이름 (length nf).
    speed_cpd : ndarray, shape (nf,)
        분조 각속도 [cycles per day].
    nodal_mode : {'amp_phase', 'phase_only', 'none'}
        nodal 보정 방식:
        'amp_phase' – f, u 모두 적용 (기본);
        'phase_only' – u만 적용 (f=1);
        'none' – nodal 보정 없음.
    R : float
        관측 잡음 분산. 결측이 많을수록 키우면 수렴이 안정적.
    P0 : float
        초기 상태 공분산 대각값. 충분히 크면 초기값 무관하게 수렴.
    store_stride : int or None
        이력 저장 간격. None이면 ~1000 스냅샷.
    epoch_dnum : float or None
        Origin of the cosine/sine arguments (MATLAB datenum). None →
        ``time_dnum[0]``. Resume from a checkpoint must pass the original
        calendar epoch, not the first time of the new slice.
    m_init, P_init : ndarray or None
        Parameter mean and covariance after the previous sample. None →
        ``m=0``, ``P=P0 I``. A checkpoint is the pair ``(m, P)``; ``m``
        alone cannot be resumed.

    Returns
    -------
    KTideResult
        result.gap_fraction : 결측 비율 (0~1)
        result.n_obs        : 실제 사용된 관측 수

    Notes
    -----
    결측 처리:
      - 관측이 NaN인 스텝: P, m 모두 업데이트하지 않음
        → 이전 추정값 유지, 불확실도 정직하게 반영
      - 보간이나 갭 메우기 불필요
      - 비등간격 관측 시계열에도 적용 가능
        (단, nodal 보정은 실제 관측 시각 기준으로 계산)
      - 이어 붙이려면 ``(m, P)`` 와 원래 ``epoch_dnum`` 을 함께 저장
    """
    time_dnum = np.asarray(time_dnum, dtype=float).ravel()
    zeta_m = np.asarray(zeta_m, dtype=float).ravel()
    speed_cpd = np.asarray(speed_cpd, dtype=float).ravel()

    nt = len(time_dnum)
    nf = len(const_names)
    t0 = float(time_dnum[0] if epoch_dnum is None else epoch_dnum)
    t_rel = time_dnum - t0   # days from the calendar epoch of G_k

    f_nc, u_nc, N_deg = compute_nodal(const_names, time_dnum)
    theta = 2.0 * np.pi * np.outer(t_rel, speed_cpd) + u_nc  # (nt, nf)

    nodal_mode = nodal_mode.lower()
    if nodal_mode == 'none':
        basis_cos = np.cos(2.0 * np.pi * np.outer(t_rel, speed_cpd))
        basis_sin = np.sin(2.0 * np.pi * np.outer(t_rel, speed_cpd))
    elif nodal_mode == 'phase_only':
        basis_cos = np.cos(theta)
        basis_sin = np.sin(theta)
    elif nodal_mode == 'amp_phase':
        basis_cos = f_nc * np.cos(theta)
        basis_sin = f_nc * np.sin(theta)
    else:
        raise ValueError(f"Unknown nodal_mode: '{nodal_mode}'")

    # ----- Kalman 초기화 -----
    if P_init is None:
        P = P0 * np.eye(2 * nf)
    else:
        P = np.array(P_init, dtype=float, copy=True)
        if P.shape != (2 * nf, 2 * nf):
            raise ValueError(f"P_init shape {P.shape} != {(2 * nf, 2 * nf)}")
    if m_init is None:
        m = np.zeros(2 * nf)
    else:
        m = np.array(m_init, dtype=float, copy=True)
        if m.shape != (2 * nf,):
            raise ValueError(f"m_init shape {m.shape} != {(2 * nf,)}")

    if store_stride is None:
        store_stride = max(1, nt // 1000)
    else:
        store_stride = max(1, int(store_stride))

    n_store = (nt - 1) // store_stride + 1 + 2
    amp_hist = np.full((n_store, nf), np.nan)
    phase_hist = np.full((n_store, nf), np.nan)
    time_hist = np.full(n_store, np.nan)
    store_idx = 0
    data_length = 0  # 격자(정점) 자료길이: 있으면 +1, 없으면 +0

    # ----- 순차 업데이트 -----
    for k in range(nt):
        obs = zeta_m[k]

        if np.isnan(obs):
            # 결측: +0, P/m 업데이트 없음
            pass
        else:
            Gk = np.empty(2 * nf)
            Gk[0::2] = basis_cos[k]
            Gk[1::2] = basis_sin[k]

            innov = obs - Gk @ m
            S = Gk @ P @ Gk + R
            K = (P @ Gk) / S
            P = P - np.outer(K, Gk @ P)
            m = m + K * innov
            data_length += 1  # 자료 있음 → +1

        if k == 0 or k == nt - 1 or k % store_stride == 0:
            a, ph = _cart_to_amp_phase(m)
            amp_hist[store_idx] = a
            phase_hist[store_idx] = ph
            time_hist[store_idx] = time_dnum[k]
            store_idx += 1

    amp_final, phase_final = _cart_to_amp_phase(m)
    n_kalman_steps = nt
    frac = data_length / n_kalman_steps if n_kalman_steps > 0 else float("nan")
    gap_fraction = 1.0 - frac if n_kalman_steps > 0 else np.nan

    return KTideResult(
        names=list(const_names),
        speed_cpd=speed_cpd.copy(),
        nodal_mode=nodal_mode,
        amp_m=amp_final,
        phase_deg=phase_final,
        model_vector=m.copy(),
        posterior_covariance=P.copy(),
        time_hist=time_hist[:store_idx],
        amp_hist=amp_hist[:store_idx],
        phase_hist_deg=phase_hist[:store_idx],
        node_longitude_deg=N_deg,
        n_obs=data_length,
        gap_fraction=gap_fraction,
        data_length=data_length,
        n_kalman_steps=n_kalman_steps,
        data_fraction=float(frac) if np.isfinite(frac) else float("nan"),
    )


# =====================================================================
# Vector (UV current) harmonic analysis
# =====================================================================

@dataclass
class KTideVectorResult:
    """Harmonic analysis results for U and V current components."""
    u: KTideResult
    v: KTideResult
    names: List[str]
    speed_cpd: np.ndarray
    nodal_mode: str
    # Tidal ellipse parameters (per constituent)
    semi_major: np.ndarray   # semi-major axis amplitude
    semi_minor: np.ndarray   # semi-minor axis amplitude (negative = CW rotation)
    inclination_deg: np.ndarray  # orientation of major axis [deg, 0-360]
    phase_deg: np.ndarray    # phase of max current along major axis [deg, 0-360]


def _compute_ellipse_params(u_amp, u_phase_deg, v_amp, v_phase_deg):
    """
    Convert U/V harmonic constants to tidal ellipse parameters.

    Uses the rotary-component method:
      w+ = (U - iV)/2  (counter-clockwise)
      w- = (U + iV)/2  (clockwise)
    where U = u_amp * exp(-i * u_phase), V = v_amp * exp(-i * v_phase).
    """
    nf = len(u_amp)
    semi_major = np.empty(nf)
    semi_minor = np.empty(nf)
    inclination = np.empty(nf)
    phase = np.empty(nf)

    for j in range(nf):
        up = np.deg2rad(u_phase_deg[j])
        vp = np.deg2rad(v_phase_deg[j])
        # complex amplitudes
        U = u_amp[j] * np.exp(-1j * up)
        V = v_amp[j] * np.exp(-1j * vp)
        # rotary components
        wp = (U + 1j * V) / 2.0   # CCW
        wm = (U - 1j * V) / 2.0   # CW (conjugate convention)
        Ap = np.abs(wp)
        Am = np.abs(wm)
        semi_major[j] = Ap + Am
        semi_minor[j] = Ap - Am   # negative → CW-dominant
        # orientation and phase
        theta_p = np.angle(wp)
        theta_m = np.angle(wm)
        inclination[j] = np.mod(np.degrees((theta_m - theta_p) / 2.0), 360.0)
        phase[j] = np.mod(np.degrees(-(theta_m + theta_p) / 2.0), 360.0)

    return semi_major, semi_minor, inclination, phase


def fit_ktide_uv(
    time_dnum: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
    const_names: List[str],
    speed_cpd: np.ndarray,
    *,
    nodal_mode: str = 'amp_phase',
    R: float = 1e-3,
    P0: float = 1e5,
    store_stride: Optional[int] = None,
) -> KTideVectorResult:
    """
    Harmonic analysis for vector (U, V) current time-series.

    U and V are analysed independently; results are combined into
    tidal ellipse parameters (semi-major, semi-minor, inclination, phase).

    Parameters
    ----------
    time_dnum : ndarray (nt,)
    u, v : ndarray (nt,)  — east/north current components [m/s or cm/s]
    const_names, speed_cpd, nodal_mode, R, P0, store_stride :
        Same as fit_ktide.

    Returns
    -------
    KTideVectorResult
    """
    # U, V 모두 NaN인 스텝은 결측 처리 (각각 독립적으로 NaN 허용)
    kw = dict(nodal_mode=nodal_mode, R=R, P0=P0, store_stride=store_stride)
    result_u = fit_ktide(time_dnum, u, const_names, speed_cpd, **kw)
    result_v = fit_ktide(time_dnum, v, const_names, speed_cpd, **kw)

    smaj, smin, inc, ph = _compute_ellipse_params(
        result_u.amp_m, result_u.phase_deg,
        result_v.amp_m, result_v.phase_deg,
    )

    return KTideVectorResult(
        u=result_u,
        v=result_v,
        names=list(const_names),
        speed_cpd=speed_cpd.copy(),
        nodal_mode=nodal_mode,
        semi_major=smaj,
        semi_minor=smin,
        inclination_deg=inc,
        phase_deg=ph,
    )


def speed_dir_to_uv(speed: np.ndarray, direction_deg: np.ndarray):
    """
    Convert current speed/direction (oceanographic convention: direction
    the current flows TOWARD, clockwise from north) to U (east) / V (north).
    """
    d = np.deg2rad(direction_deg)
    u = speed * np.sin(d)
    v = speed * np.cos(d)
    return u, v
