"""
nodal_corrections.py
Approximate nodal amplitude (f) and phase (u) corrections.
Ported from khoa_compute_nodal_extended.m
"""
import numpy as np

# ---------------------------------------------------------------------------
# Group mapping: constituent name → base group for nodal table look-up
# ---------------------------------------------------------------------------
_M2_GROUP = {
    'M2','N2','2N2','MU2','NU2','L2','LAM2','MKS2','MSN2','OQ2',
    'MNS2','OP2','MA2','MB2',
}
_K2_GROUP = {'K2', 'KJ2'}
_K1_GROUP = {'K1', 'P1', 'PI1', 'PSI1', 'PHI1', 'S1'}
_O1_GROUP = {
    'O1','Q1','2Q1','SIG1','RO1','MP1','M1','CHI1','TH1','J1','SO1','OO1',
}
_MF_GROUP = {'MF', 'MSF'}
_MM_GROUP = {'MM'}

# base nodal table: {name: (f_coef, u_coef)}
# f(N) = f_coef[0] + f_coef[1]*cos(N) + f_coef[2]*cos(2N) + ...
# u(N) = u_coef[0]*sin(N) + u_coef[1]*sin(2N) + ...  [degrees]
_BASE_TABLE = {
    'M2': ([1.0004, -0.0373,  0.0002], [-2.14,  0.10]),
    'S2': ([1.0],                       [0.0]),
    'K2': ([1.0241,  0.2863,  0.0083], [-17.74, 0.68]),
    'K1': ([1.0060,  0.1150, -0.0088,  0.0006], [-8.86, 0.68, -0.07]),
    'O1': ([1.0089,  0.1871, -0.0147], [10.80, -1.34]),
    'MF': ([1.0429,  0.4135, -0.0040], [-23.74, 2.68]),
    'MM': ([1.0000, -0.1300],           [0.0]),
}


def _group_name(name: str) -> str:
    n = name.upper().strip()
    if n in _M2_GROUP:
        return 'M2'
    if n in _K2_GROUP:
        return 'K2'
    if n in _K1_GROUP:
        return 'K1'
    if n in _O1_GROUP:
        return 'O1'
    if n in _MF_GROUP:
        return 'MF'
    if n in _MM_GROUP:
        return 'MM'
    return n


def _lunar_node_deg(time_dnum: np.ndarray) -> np.ndarray:
    """Lunar ascending node longitude [deg] from MATLAB datenums."""
    JD = time_dnum + 1721058.5
    T = (JD - 2451545.0) / 36525.0
    N = (125.0445479
         - 1934.1362608 * T
         + 0.0020754 * T**2
         + T**3 / 467441.0
         - T**4 / 60616000.0)
    return N % 360.0


def compute_nodal(const_names: list, time_dnum: np.ndarray):
    """
    Compute nodal amplitude factors (f) and phase offsets (u) for each
    constituent at each time step.

    Parameters
    ----------
    const_names : list of str
        Constituent names, length nf.
    time_dnum : ndarray, shape (nt,)
        MATLAB serial datenums.

    Returns
    -------
    f_nc : ndarray, shape (nt, nf)
        Nodal amplitude factors (dimensionless, ~1).
    u_nc : ndarray, shape (nt, nf)
        Nodal phase offsets [radians].
    N_deg : ndarray, shape (nt,)
        Lunar node longitude [degrees].
    """
    time_dnum = np.asarray(time_dnum, dtype=float).ravel()
    nt = len(time_dnum)
    nf = len(const_names)

    N_deg = _lunar_node_deg(time_dnum)
    N_rad = np.deg2rad(N_deg)

    f_nc = np.ones((nt, nf))
    u_nc = np.zeros((nt, nf))

    for j, name in enumerate(const_names):
        base = _group_name(name)
        if base not in _BASE_TABLE:
            continue
        fc, uc = _BASE_TABLE[base]

        fj = fc[0] * np.ones(nt)
        for n, c in enumerate(fc[1:], start=1):
            fj += c * np.cos(n * N_rad)

        uj = np.zeros(nt)
        for n, c in enumerate(uc, start=1):
            uj += c * np.sin(n * N_rad)

        f_nc[:, j] = fj
        u_nc[:, j] = np.deg2rad(uj)

    return f_nc, u_nc, N_deg
