# ─────────────────────────────────────────────────────────────────
# Cas13 gRNA Designer — backend Dockerfile (build context = the backend folder)
# AI: Groq API (no GPU needed)
# Includes: ViennaRNA, RNAstructure, LinearFold, Python, MySQL client
# ─────────────────────────────────────────────────────────────────

FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

# ─── SYSTEM PACKAGES ──────────────────────────────────────────────
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-dev \
    build-essential \
    wget \
    curl \
    git \
    cmake \
    libssl-dev \
    libffi-dev \
    pkg-config \
    default-mysql-client \
    default-libmysqlclient-dev \
    && rm -rf /var/lib/apt/lists/*

# ─── INSTALL ViennaRNA ────────────────────────────────────────────
RUN apt-get update && apt-get install -y \
    vienna-rna \
    && rm -rf /var/lib/apt/lists/*

RUN RNAfold --version && echo "RNAfold OK"
# RNAplfold gives seed confidence for long sequences; fail the build early if it is missing
RUN RNAplfold --version && echo "RNAplfold OK"

# ─── INSTALL RNAstructure ─────────────────────────────────────────
# Optional tool: if the download fails the build continues and the app
# simply reports RNAstructure as unavailable (RNAfold still works).
RUN (wget -q \
    "https://rna.urmc.rochester.edu/Releases/current/RNAstructureLinuxTextInterfaces.tgz" \
    -O /tmp/rnastructure.tgz && \
    tar -xzf /tmp/rnastructure.tgz -C /opt/ && \
    rm /tmp/rnastructure.tgz) || echo "WARNING: RNAstructure install failed - continuing without it"

ENV PATH="/opt/RNAstructure/exe:${PATH}"
ENV DATAPATH="/opt/RNAstructure/data_tables"

RUN Fold --version || echo "RNAstructure installed"

# ─── INSTALL LinearFold ───────────────────────────────────────────
# Optional tool: the app falls back to windowed RNAfold if this is missing.
# We put the repo folder on PATH instead of copying one file out of it,
# because the 'linearfold' launcher looks for its compiled binaries in bin/.
RUN (git clone --depth=1 \
    https://github.com/LinearFold/LinearFold.git \
    /opt/linearfold && \
    cd /opt/linearfold && \
    make && \
    chmod +x /opt/linearfold/linearfold && \
    echo "LinearFold installed") || echo "WARNING: LinearFold build failed - app will use windowed fallback"

ENV PATH="/opt/linearfold:${PATH}"

# ─── PYTHON APP ───────────────────────────────────────────────────
WORKDIR /app

COPY requirements.txt .
RUN pip3 install --no-cache-dir -r requirements.txt

COPY . .

# ─── ENV DEFAULTS ─────────────────────────────────────────────────
ENV FLASK_DEBUG=False
ENV FLASK_PORT=5000
ENV DB_HOST=db
ENV DB_USER=root
ENV DB_NAME=cas13_grna_db
ENV MAX_SEQUENCE_LENGTH=20000
# Groq model: can be overridden by a Railway variable of the same name
ENV GROQ_MODEL=openai/gpt-oss-20b

EXPOSE 5000

COPY start.sh /start.sh
# strip Windows CRLF line endings (they make the script fail with "no such file")
RUN sed -i 's/\r$//' /start.sh && chmod +x /start.sh

CMD ["/start.sh"]
