# GitHub + Zenodo (DOI)

Do this after `/pao1/work/papers/ktide/deposit` is the GitHub root.
This folder is separate from `papers/cycloEOF/deposit/` (nested CSEOF).
Do not mix the two archives. MATLAB KTide, T_TIDE, UTide, TIRA, and
Incheon gauge records are not in this folder. Do not add them.

## 0. GitHub repo (account syongkim)

On the machine with git and SSH to GitHub:

```bash
cd /pao1/work/papers/ktide/deposit
git init
git add LICENSE NOTICE README.md CITATION.cff ZENODO.md requirements.txt \
  .gitignore ktide examples tables
git commit -m "KTide Python: Q=0 sequential harmonic analysis as distributed."
```

Create an empty public repository `KTide` under https://github.com/syongkim
(no README on GitHub, so the first push is clean). Then:

```bash
git branch -M main
git remote add origin git@github.com:syongkim/KTide.git
git push -u origin main
```

Zenodo's GitHub harvest needs a **public** repo (or a paid Zenodo/GitHub
connection). Keep the repo public if the paper's availability statement
should resolve without a login.

## 1. Zenodo account

1. Open https://zenodo.org and **Log in with GitHub** (same `syongkim`).
   That is the least error-prone link. ORCID login also works; then
   connect GitHub under the account menu.
2. Confirm the email Zenodo sends.
3. You do not need to share a password. After login, the GitHub
   applications list will show Zenodo.

## 2. Flip the repo on

1. Zenodo: GitHub icon (top right) → **GitHub**.
2. Find `syongkim/KTide` and switch it **on**.
3. Zenodo installs a webhook on that repository.

If the repo does not appear, grant Zenodo access to that repository in
GitHub → Settings → Applications → Zenodo → Repository access.

## 3. Mint the DOI (a GitHub Release)

Zenodo does not mint a DOI from an ordinary push. It mints one when
GitHub publishes a **Release**.

1. GitHub → `KTide` → **Releases** → **Draft a new release**.
2. Tag: `v1.0.0` (create tag on `main`).
3. Title: `KTide v1.0.0`.
4. Publish release.

Wait one to several minutes. Zenodo harvests the tag, builds a tarball,
and assigns:

- a **version DOI** for `v1.0.0` (cite this in the paper), e.g.
  `10.5281/zenodo.1234567`
- a **concept DOI** that always points at the latest version (use in
  README badges).

Record both. Put the version DOI in `CITATION.cff` (`doi:`) and in
`\datastatement` of `ktide_jtech.tex`. Later tagged releases get new
version DOIs; the concept DOI stays.

## 4. What not to upload

- `ksstats` (separate repository `syongkim/ksstats`)
- `papers/ktide/python/` experiment drivers and `out/`
- MATLAB KTide (`ktide_setup.m` and the operational tree)
- T_TIDE / UTide / TIRA distributions
- Incheon hourly gauge records
- `papers/cycloEOF/deposit/` (different manuscript)

If a file is not needed to run `examples/` and to read the comparison
CSVs in `tables/`, leave it out.
