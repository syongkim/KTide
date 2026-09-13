"""
constituents.py
Standard constituent set used in Korea (KHOA) tidal analysis.
Ported from khoa_tide_constituents_case_incheon.m
"""
import numpy as np
from dataclasses import dataclass, field
from typing import List


@dataclass
class ConstituentSet:
    names: List[str]
    speed_cph: np.ndarray   # deg/hour
    speed_cpd: np.ndarray   # cycles/day  (= speed_cph / 15)
    amp_m: np.ndarray       # amplitude [m]   (benchmark truth values)
    phase_deg: np.ndarray   # phase [deg, 0-360]
    source: str = ""


def incheon_constituents() -> ConstituentSet:
    """
    Representative synthetic constituent set for a Korea west-coast case
    (Incheon).  Amplitudes and phases are seeded from published harmonic
    constants and are intended for benchmark use, not operational prediction.

    Ported from khoa_tide_constituents_case_incheon.m
    """
    names = [
        'Z0','SA','SSA','MM','MSF','MF',
        '2Q1','SIG1','Q1','CHI1','PI1','P1','S1','K1',
        'PSI1','PHI1','J1','SO1','OO1',
        'OQ2','2N2','MU2','N2','NU2','OP2','M2','MKS2','L2',
        'T2','S2','R2','K2','MSN2','2SM2',
        'MO3','M3','SO3','MK3','SK3',
        'MN4','M4','SN4','MS4','MK4','S4','SK4',
        '2MN6','M6','MSN6','2MS6','2MK6','2SM6','MSK6',
    ]

    speed_cph = np.array([
        0.0000000, 0.0410686, 0.0821373, 0.5443747, 1.0158958, 1.0980331,
        12.8542862, 12.9271398, 13.3986609, 14.5695476, 14.9178647, 14.9589314,
        15.0000000, 15.0410686, 15.0821353, 15.1232059,
        15.5854433, 16.0569644, 16.1391017,
        27.3416965, 27.8953548, 27.9682084, 28.4397295,
        28.5125831, 28.9019669, 28.9841042, 29.0662415, 29.5284789,
        29.9589333, 30.0000000, 30.0410667, 30.0821373,
        30.5443747, 31.0158958,
        42.9271398, 43.4761563, 43.9430356, 44.0251729, 45.0410686,
        57.4238337, 57.9682084, 58.4397295, 58.9841042, 59.0662415,
        60.0000000, 60.0821373,
        86.4079380, 86.9523127, 87.4238337, 87.9682084, 88.0503457,
        88.9841042, 89.0662415,
    ])

    amp_m = np.array([
        0.000, 0.080, 0.040, 0.030, 0.060, 0.090,
        0.006, 0.008, 0.047, 0.008,
        0.007, 0.129, 0.010, 0.389, 0.003, 0.005, 0.020, 0.008, 0.012,
        0.006, 0.070, 0.094, 0.525, 0.102, 0.012, 2.823, 0.015,
        0.150, 0.067, 1.135, 0.030, 0.309, 0.012, 0.025,
        0.006, 0.040, 0.008, 0.012, 0.005,
        0.025, 0.100, 0.020, 0.070, 0.015, 0.030, 0.010,
        0.008, 0.040, 0.015, 0.012, 0.010, 0.008, 0.008,
    ])

    phase_deg = np.mod(np.array([
          0.0, 200.0, 230.0,  20.0,  35.0,  40.0,
        230.0, 246.0, 238.7, 288.0,
        296.0, 300.0, 301.0, 303.2, 306.0, 309.0, 322.0, 334.0, 340.0,
         96.0, 104.0, 178.3, 108.4, 112.0, 120.0, 129.5, 138.0,
        148.8, 180.0, 186.1, 190.0, 192.0, 205.0, 220.0,
         35.0,  42.0,  55.0,  60.0,  72.0,
        258.0, 266.0, 280.0, 290.0, 302.0, 315.0, 325.0,
         28.0,  38.0,  48.0,  58.0,  68.0,  78.0,  88.0,
    ]), 360.0)

    return ConstituentSet(
        names=names,
        speed_cph=speed_cph,
        speed_cpd=speed_cph / 15.0,
        amp_m=amp_m,
        phase_deg=phase_deg,
        source=(
            'Synthetic benchmark seeded from representative '
            'Incheon constituent magnitudes'
        ),
    )
