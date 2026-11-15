# Backend (Flask API, deployed on Railway)

| File | Purpose |
|---|---|
| `app_py_final.py` | The API: `/` health, `/analyze`, `/ai_agent`, `/ai_chat` |
| `linearfold_integration.py` | Structure prediction (RNAfold, RNAstructure, LinearFold, RNAplfold) |
| `setup_db.py` | Creates the MySQL tables (`--reset` drops and recreates them) |
| `load_pars_real.py` | Loads the experimental PARS scores into MySQL |
| `load_literature.py` | Loads the guides and single-stranded regions published in Bandaru et al. 2020 |
| `start.sh` | Container start-up: wait for MySQL, set up, load data, start gunicorn |
| `Dockerfile` | Installs ViennaRNA, LinearFold, optional RNAstructure, Python libraries |
| `pars_real_GM12878.tsv.gz` | The real PARS data. **This copy holds XIST only (a starter).** Replace it with your full file (62 transcripts) before deploying |
| `literature/` | Bandaru et al. 2020 guide list and single-stranded regions for XIST |
| `tools/build_pars_database.py` | Rebuilds `pars_real_GM12878.tsv.gz` from the GEO files |
| `benchmark_bandaru2020.py` | Re-run the cross-check against the paper after any change |

## Environment variables (Railway: Variables tab)

`DB_HOST`, `DB_PORT`, `DB_USER`, `DB_PASSWORD`, `DB_NAME` (set to `cas13_grna_db`), `GROQ_API_KEY`, `GROQ_MODEL`.
See `.env.example`. `PORT` is set by Railway automatically.

## Run locally

```
pip install -r requirements.txt
python setup_db.py --reset
python load_pars_real.py
python load_literature.py
python app_py_final.py
```
