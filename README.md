# PISA 2025 generative-AI frequency paradox — analysis code

This repository contains the analysis scripts, variable crosswalk and download log for the manuscript
"The Frequency Paradox of Generative AI: Task-Contingent Curvature and AI-Literacy Moderation in Science Achievement".

## Data

The analysis uses the OECD PISA 2025 Student Public Use File (PUF), which is **not redistributed here**
because of its size and licence. Download it from the official OECD PISA data site; the exact file name,
URL and integrity checks are recorded in `DOWNLOAD_LOG.md`.

Required official files:

- `CY09_MS_STU_PUF.sav` (student microdata)
- PISA 2025 Codebook and student/ICT questionnaires (for variable mapping)

## Contents

| File | Purpose |
|---|---|
| `run_real_analysis.py` | Main WLS / BRR / Rubin pooled estimates |
| `G04_forensics.py` | Component/forensic checks |
| `formal_by_country.py` | Country-level estimates and random-effects meta-analysis |
| `meta_regression.py` | Random-effects meta-regression moderators |
| `supplementary_analyses.py` | H1a task-specific, H2b falsification, H3/H4 stratification |
| `h2b_joint_test.py` | Joint AI-specific vs general-online contrast |
| `method_formal_tests.py` | Lind–Mehlum, functional forms, Johnson–Neyman, FDR |
| `generate_manuscript_figures.py` | Produces figures 2–12 (SVG + PNG) |
| `final_reaudit.py` | Read-only consistency and disclosure audit of the manuscript |
| `codebook_ai_pisa2025.csv` | Variable crosswalk used for item mapping |
| `public_indicators_from_sources.csv` | Official public indicators transcribed from source tables |
| `DOWNLOAD_LOG.md` | Official download URL, file integrity and version checks |

## Requirements

Python 3.10+ with `numpy`, `scipy`, `pandas`, `matplotlib`, `pyreadstat`.

## Reproducing the reported estimates

1. Download the official PISA 2025 student PUF and supporting documents as described in `DOWNLOAD_LOG.md`.
2. Place the `.sav` file in the path expected by `run_real_analysis.py` and run the scripts in this order:
   `run_real_analysis.py`, `G04_forensics.py`, `formal_by_country.py`, `meta_regression.py`,
   `supplementary_analyses.py`, `h2b_joint_test.py`, `method_formal_tests.py`,
   `generate_manuscript_figures.py`, `final_reaudit.py`.
3. `final_reaudit.py` checks manuscript consistency and disclosure; it does **not** refit models or alter results.

All estimates in the manuscript are cross-sectional associations. Vertices are turning points of fitted
cross-sectional curves and are not intervention optima.

## Licence

MIT. See `LICENSE`.
