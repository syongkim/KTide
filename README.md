# KTide (Python)

Sequential recursive least-squares harmonic analysis (\(Q=0\)) of
stored coastal archives, as recorded in Kim and Byun (in prep.,
*J. Atmos. Oceanic Technol.*).

This is the paper implementation. MATLAB KTide is **not** included.
T\_TIDE, UTide, and TIRA are **not** included. Sequential moments
(`KSMoments`, `KSPairs`) live in a **separate** repository,
[syongkim/ksstats](https://github.com/syongkim/ksstats). See `NOTICE`.

## Run

```bash
python -m pip install -r requirements.txt
PYTHONPATH=. python examples/fit_one_series.py
PYTHONPATH=. python examples/resume_checkpoint.py
```

`fit_ktide` times are MATLAB serial datenum (days). The example builds
them from `numpy.datetime64`. Resume from a checkpoint by passing the
original calendar epoch (`epoch_dnum`) together with `(m_init, P_init)`;
`m` alone cannot be resumed.

## Contents

| Path | What |
|------|------|
| `ktide/` | Sequential recursive least squares, \(Q=0\) (`fit_ktide`, `fit_ktide_uv`) |
| `examples/` | Minimal calls matching Appendix B of the manuscript |
| `tables/` | Comparison CSVs used for the four-method figures and Table A1; `tables/v1.1.0/` adds the tables of the archive, unstructured-mesh, D1 ensemble, uncertainty, and file-read experiments |

Not included: `ksstats`, four-method driver runs, Incheon hourly gauge
records, KHOA operational wrappers.

## License

MIT (`LICENSE`) for this Python. T\_TIDE and UTide remain under their
authors' terms. Method literature remains with its authors.

## Cite

Version DOI of `v1.1.0`: [10.5281/zenodo.23051866](https://doi.org/10.5281/zenodo.23051866).
Version DOI of `v1.0.0`: [10.5281/zenodo.22738041](https://doi.org/10.5281/zenodo.22738041).

## v1.1.0 (2026-09-30)

- `ktide.fit_ktide`: optional increment variance `Q` (default 0, recursive least squares).
- `ktide.fit_ktide_shared`: shared-covariance map-slice update in double precision for cells sampled at the same times; refuses missing samples (cells with different masks need their own `P`).
- `tables/v1.1.0/`: numbers behind the stored-archive (YES3k, MOHID), SCHISM wetting-drying, D1 Monte Carlo (1000 realizations), uncertainty, and file-read experiments.
- Title follows the revised paper (no "Sequential").
The concept DOI [10.5281/zenodo.22738040](https://doi.org/10.5281/zenodo.22738040)
always points at the latest version. Cite the paper for the method.
