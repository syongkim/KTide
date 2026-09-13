"""
benchmark.py
Build a synthetic tidal time-series and compare KTide (and optionally
UTide / T_TIDE) harmonic analysis results against the prescribed truth.
Ported from khoa_run_tide_benchmark.m
"""
from __future__ import annotations

import os
import warnings
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd

from .constituents import ConstituentSet, incheon_constituents
from .nodal_corrections import compute_nodal
from .fit import fit_ktide, KTideResult


# ---------------------------------------------------------------------------
# Helper: MATLAB datenum ↔ datetime
# ---------------------------------------------------------------------------
_MATLAB_EPOCH = np.datetime64('0000-01-01') - np.timedelta64(1, 'D')

def datetime64_to_dnum(dt64: np.ndarray) -> np.ndarray:
    """Convert numpy datetime64[s] → MATLAB datenum (float64)."""
    delta = (dt64 - _MATLAB_EPOCH).astype('timedelta64[s]').astype(float)
    return delta / 86400.0

def dnum_to_datetime64(dnum: np.ndarray) -> np.ndarray:
    """Convert MATLAB datenum → numpy datetime64[s]."""
    sec = np.round((dnum) * 86400).astype('int64')
    return _MATLAB_EPOCH + sec.astype('timedelta64[s]')


# ---------------------------------------------------------------------------
# Build synthetic truth series
# ---------------------------------------------------------------------------
def build_truth(
    C: ConstituentSet,
    *,
    start: str = '2000-01-01',
    years: float = 22.0,
    dt_hours: float = 1.0,
    noise_std_m: float = 0.0,
    seed: int = 42,
) -> dict:
    """
    Generate a synthetic tidal time-series from prescribed constituents.

    Returns a dict with keys:
      time_dnum, zeta_truth_m, zeta_m, amp_m, phase_deg, names, speed_cpd
    """
    t0 = np.datetime64(start, 's')
    n_steps = int(round(years * 365.2425 * 24 / dt_hours)) + 1
    dt_sec = int(dt_hours * 3600)
    time_dt = t0 + np.arange(n_steps, dtype='int64') * np.timedelta64(dt_sec, 's')
    time_dnum = datetime64_to_dnum(time_dt)

    f_nc, u_nc, _ = compute_nodal(C.names, time_dnum)
    t_rel = time_dnum - time_dnum[0]
    theta = (2.0 * np.pi * np.outer(t_rel, C.speed_cpd)
             + u_nc
             - np.deg2rad(C.phase_deg))
    zeta = np.sum(f_nc * C.amp_m * np.cos(theta), axis=1)

    if noise_std_m > 0:
        rng = np.random.default_rng(seed)
        zeta = zeta + noise_std_m * rng.standard_normal(zeta.shape)

    return dict(
        time_dnum=time_dnum,
        time_dt=time_dt,
        zeta_truth_m=zeta.copy(),
        zeta_m=zeta.copy(),
        amp_m=C.amp_m.copy(),
        phase_deg=C.phase_deg.copy(),
        names=list(C.names),
        speed_cpd=C.speed_cpd.copy(),
    )


# ---------------------------------------------------------------------------
# Error table helpers
# ---------------------------------------------------------------------------
def _wrap180(x: np.ndarray) -> np.ndarray:
    return (x + 180.0) % 360.0 - 180.0


def compare_to_truth(
    C: ConstituentSet,
    est_names: list,
    est_amp: np.ndarray,
    est_phase: np.ndarray,
    method_name: str,
) -> pd.DataFrame:
    """Build a comparison table between estimated and truth constituents."""
    truth_names = [n.upper() for n in C.names]
    est_map = {n.upper(): (a, p)
               for n, a, p in zip(est_names, est_amp, est_phase)}

    rows = []
    for name, t_amp, t_phase, speed in zip(
            truth_names, C.amp_m, C.phase_deg, C.speed_cpd):
        if name in est_map:
            e_amp, e_phase = est_map[name]
            e_phase = e_phase % 360.0
        else:
            e_amp, e_phase = np.nan, np.nan
        amp_err = e_amp - t_amp if not np.isnan(e_amp) else np.nan
        amp_pct = 100.0 * amp_err / max(t_amp, 1e-6) if not np.isnan(amp_err) else np.nan
        ph_err = _wrap180(np.array([e_phase - t_phase]))[0] if not np.isnan(e_phase) else np.nan
        rows.append(dict(
            method=method_name,
            name=name,
            speed_cpd=speed,
            truth_amp_m=t_amp,
            truth_phase_deg=t_phase,
            est_amp_m=e_amp,
            est_phase_deg=e_phase,
            amp_error_m=amp_err,
            amp_error_pct=amp_pct,
            phase_error_deg=ph_err,
        ))
    return pd.DataFrame(rows)


def summarize_error(df: pd.DataFrame, method_name: str) -> dict:
    valid = df['est_amp_m'].notna()
    return dict(
        method=method_name,
        n_found=int(valid.sum()),
        n_total=len(df),
        amp_rmse_m=float(np.sqrt(np.mean(df.loc[valid, 'amp_error_m']**2))),
        phase_rmse_deg=float(np.sqrt(np.mean(df.loc[valid, 'phase_error_deg']**2))),
        amp_mean_abs_pct=float(np.mean(np.abs(df.loc[valid, 'amp_error_pct']))),
    )


# ---------------------------------------------------------------------------
# Main benchmark runner
# ---------------------------------------------------------------------------
def run_benchmark(
    *,
    years: float = 22.0,
    dt_hours: float = 1.0,
    start: str = '2000-01-01',
    noise_std_m: float = 0.0,
    results_dir: str = 'benchmark_out',
    utide_available: bool = False,
) -> dict:
    """
    Build synthetic tidal series and run KTide harmonic analysis.

    Parameters
    ----------
    years : float
        Length of the synthetic series (default 22 yr → covers nodal cycle).
    dt_hours : float
        Time step in hours.
    start : str
        Start date string (ISO format).
    noise_std_m : float
        Gaussian observation noise standard deviation [m].
    results_dir : str
        Directory to write CSV / PNG output.
    utide_available : bool
        If True and the `utide` Python package is installed, reconstruct the
        benchmark series via UTide before fitting KTide.

    Returns
    -------
    dict with keys: truth, KTide, tables, summaries
    """
    os.makedirs(results_dir, exist_ok=True)

    C = incheon_constituents()
    truth = build_truth(C, start=start, years=years,
                        dt_hours=dt_hours, noise_std_m=noise_std_m)

    zeta_bench = truth['zeta_m']
    time_bench = truth['time_dnum']

    # Optionally replace benchmark series with UTide reconstruction
    if utide_available:
        try:
            import utide
            coef = utide.solve(
                truth['time_dt'].astype('datetime64[s]').astype(float) / 86400.0,
                zeta_bench,
                lat=37.45,
                constit=C.names,
            )
            recon = utide.reconstruct(
                truth['time_dt'].astype('datetime64[s]').astype(float) / 86400.0,
                coef,
            )
            zeta_bench = recon.h
        except Exception as exc:
            warnings.warn(f"UTide reconstruction failed: {exc}. Using direct series.")

    print(f"시계열 길이: {len(time_bench):,} 스텝  ({years:.1f}년, {dt_hours}h 간격)")
    print("KTide 조화분해 시작 (nodal amp+phase only)…")

    r = fit_ktide(time_bench, zeta_bench, C.names, C.speed_cpd,
                  nodal_mode='amp_phase')
    results = {'KTide': r}
    print("완료")

    # Build comparison tables
    tables = {}
    summaries = {}
    for key, rr in results.items():
        df = compare_to_truth(C, rr.names, rr.amp_m, rr.phase_deg, key)
        tables[key] = df
        summaries[key] = summarize_error(df, key)
        df.to_csv(os.path.join(results_dir, f'compare_{key}.csv'), index=False)

    # Summary print
    print("\n=== 분조 복원 오류 요약 ===")
    for key, s in summaries.items():
        print(f"  {s['method']:30s}  "
              f"amp_RMSE={s['amp_rmse_m']:.4f} m  "
              f"phase_RMSE={s['phase_rmse_deg']:.2f}°  "
              f"({s['n_found']}/{s['n_total']} 분조)")

    # Plot
    try:
        _make_plot(C, results, results_dir)
    except Exception as exc:
        warnings.warn(f"Plot 생성 실패: {exc}")

    return dict(
        truth=truth,
        tables=tables,
        summaries=summaries,
        **results,
    )


def _make_plot(C: ConstituentSet, results: dict, results_dir: str):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.dates import AutoDateLocator, DateFormatter

    major = [n in {'Q1','O1','K1','N2','M2','S2','K2','L2','M4','MS4'}
             for n in C.names]
    labels = [n for n, m in zip(C.names, major) if m]
    truth_amp = C.amp_m[major]
    r_ap = results['KTide']

    fig, axes = plt.subplots(2, 1, figsize=(14, 8))

    # — Amplitude bar chart —
    bar_data = np.column_stack([
        truth_amp,
        r_ap.amp_m[major],
    ])
    x = np.arange(len(labels))
    w = 0.35
    for i, lbl in enumerate(['Truth', 'KTide']):
        axes[0].bar(x + (i - 0.5) * w, bar_data[:, i], w, label=lbl)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels)
    axes[0].set_ylabel('Amplitude (m)')
    axes[0].set_title('Major constituent amplitude recovery')
    axes[0].legend(ncol=2)

    # — Convergence plot —
    time_dt = dnum_to_datetime64(r_ap.time_hist).astype('datetime64[ms]')
    import matplotlib.dates as mdates
    t_plot = mdates.date2num(time_dt.astype('M8[ms]').astype('O'))
    for j, (name, flag) in enumerate(zip(C.names, major)):
        if flag:
            axes[1].plot(t_plot, r_ap.amp_hist[:, j], lw=0.8, label=name)
    axes[1].xaxis.set_major_locator(AutoDateLocator())
    axes[1].xaxis.set_major_formatter(DateFormatter('%Y'))
    axes[1].set_ylabel('Amplitude estimate (m)')
    axes[1].set_title('KTide amp+phase – convergence of major constituents')
    axes[1].legend(ncol=5, fontsize=7)

    fig.tight_layout()
    out_path = os.path.join(results_dir, 'ktide_summary.png')
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  그림 저장: {out_path}")
