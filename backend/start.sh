#!/bin/bash
# start.sh — Container startup script
# Cas13 gRNA Designer v2.0
# Waits for MySQL, sets up tables, loads data, then runs gunicorn

set -e

echo "============================================================"
echo "  Cas13 gRNA Designer v2.0 — Container Startup"
echo "============================================================"

# ─── Check required variables ─────────────────────────────────────
for v in DB_HOST DB_PASSWORD; do
    if [ -z "${!v}" ]; then
        echo "  ✗ Required environment variable $v is not set. Exiting."
        exit 1
    fi
done
if [ -z "$GROQ_API_KEY" ]; then
    echo "  ! GROQ_API_KEY is not set - AI analysis will be unavailable"
fi

# ─── Wait for MySQL (uses Python so port and auth plugin just work) ─
echo "[1/4] Waiting for MySQL to be ready..."
python3 - <<'PY'
import os, sys, time
import mysql.connector
cfg = dict(host=os.environ["DB_HOST"],
           user=os.environ.get("DB_USER", "root"),
           password=os.environ["DB_PASSWORD"],
           port=int(os.environ.get("DB_PORT", "3306")),
           connection_timeout=5)
for attempt in range(1, 61):
    try:
        mysql.connector.connect(**cfg).close()
        print("  ✓ MySQL is ready")
        sys.exit(0)
    except Exception as e:
        print("  ... waiting for MySQL (attempt %d/60): %s" % (attempt, e))
        time.sleep(2)
print("  ✗ MySQL did not become ready. Exiting.")
sys.exit(1)
PY

# ─── Setup database ───────────────────────────────────────────────
echo "[2/4] Setting up database..."
python3 setup_db.py
echo "  ✓ Database tables ready"

# ─── Load real experimental PARS data ────────────────────────────
# pars_real_GM12878.tsv.gz holds real PARS scores (GEO GSE50676). Nothing is simulated.
echo "[3/4] Loading real PARS data into the database..."
python3 load_pars_real.py
echo "  ✓ PARS data loaded"
echo "      Loading published guides from Bandaru et al. 2020 (XIST)..."
python3 load_literature.py

# ─── Start server ─────────────────────────────────────────────────
echo "[4/4] Starting gunicorn on port ${PORT:-5000}..."
echo "============================================================"

exec gunicorn \
    --bind 0.0.0.0:${PORT:-5000} \
    --workers 2 \
    --timeout 600 \
    --keep-alive 5 \
    --log-level info \
    --access-logfile - \
    --error-logfile - \
    app_py_final:app
