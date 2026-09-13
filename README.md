# KTide (Python)

Sequential Kalman harmonic analysis (\(Q=0\) recursive least squares) of
stored coastal archives, as recorded in Kim and Byun (submitted to
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
| `ktide/` | Sequential Kalman harmonic analysis (`fit_ktide`, `fit_ktide_uv`) |
| `examples/` | Minimal calls matching Appendix B of the manuscript |
| `tables/` | Comparison CSVs used for the ranked figures and Table A1 |

Not included: `ksstats`, four-method driver runs, Incheon hourly gauge
records, KHOA operational wrappers.

## License

MIT (`LICENSE`) for this Python. T\_TIDE and UTide remain under their
authors' terms. Method literature remains with its authors.

## Cite

Version DOI of `v1.0.0`: [10.5281/zenodo.22738041](https://doi.org/10.5281/zenodo.22738041).
The concept DOI [10.5281/zenodo.22738040](https://doi.org/10.5281/zenodo.22738040)
always points at the latest version. Cite the paper for the method.
