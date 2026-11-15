# Cas13 gRNA Designer

Design Cas13 guide RNAs that target single-stranded RNA, following Bandaru et al. 2020 (Sci Rep 10:11610).
Guides are checked against real PARS experiments (GEO GSE50676) where data exists.

```
backend/    Flask API + data loaders  -> deploy on Railway (Root Directory = backend)
frontend/   Static website            -> deploy on Netlify (Publish directory = frontend)
```

## Deploy

1. **Backend (Railway).** New Project > Deploy from GitHub repo. In the service Settings set
   **Root Directory = `backend`** (Railway then finds the Dockerfile). Add a MySQL database to the project.
   Variables: `DB_HOST=${{MySQL.MYSQLHOST}}`, `DB_PORT=${{MySQL.MYSQLPORT}}`, `DB_USER=${{MySQL.MYSQLUSER}}`,
   `DB_PASSWORD=${{MySQL.MYSQLPASSWORD}}`, `DB_NAME=cas13_grna_db`, `GROQ_API_KEY=...`, `GROQ_MODEL=openai/gpt-oss-20b`.
   Settings > Networking > Generate Domain.
2. **Frontend (Netlify).** Edit `frontend/index.html`: set `window.CAS13_API` to the Railway URL (https, no trailing slash),
   commit, then import the repo on Netlify with **no build command** and **Publish directory = `frontend`**.
3. **Check.** Open the Railway URL: it should show `"database": "connected"` and the number of transcripts loaded.

Predicted guides are a shortlist: structure is predicted without experimental constraints and agrees only weakly with
lab results on XIST. Test several top guides at the bench.
