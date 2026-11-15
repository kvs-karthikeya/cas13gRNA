# Cas13 gRNA Designer — Railway backend
# Build context: backend/
# Includes: Python, MySQL client, ViennaRNA, RNAstructure, LinearFold

FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

# SYSTEM DEPENDENCIES
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
    libgomp1 \
    default-mysql-client \
    default-libmysqlclient-dev \
    vienna-rna \
    && rm -rf /var/lib/apt/lists/*

# VIENNARNA (required: RNAfold folds, RNAplfold gives seed confidence for long sequences)
RUN RNAfold --version && echo "RNAfold OK"
RUN RNAplfold --version && echo "RNAplfold OK"

# RNASTRUCTURE 64-BIT TEXT INTERFACES (optional)
# Tries several official mirrors. If none works the build still succeeds with a clear warning
# and the app uses RNAfold instead; a temporary outage of the Rochester server must not stop a deploy.
ARG RNASTRUCTURE_VERSION=6.5
RUN set -e; \
    for host in rna.urmc.rochester.edu rna2.urmc.rochester.edu rna1.urmc.rochester.edu; do \
      echo "Trying ${host}..."; \
      if timeout 45s wget --timeout=15 --tries=1 \
          "https://${host}/Releases/6.5/RNAstructureLinuxTextInterfaces64bit.tgz" \
          -O /tmp/rnastructure.tgz \
          && tar -xzf /tmp/rnastructure.tgz -C /opt/ \
          && test -x /opt/RNAstructure/exe/Fold; then \
        echo "RNAstructure 6.5 installed from ${host}"; \
        break; \
      fi; \
      rm -f /tmp/rnastructure.tgz; \
    done; \
    if test -x /opt/RNAstructure/exe/Fold; then \
      echo "RNAstructure available"; \
    else \
      echo "RNAstructure unavailable; using ViennaRNA RNAfold fallback"; \
    fi

ENV PATH="/opt/RNAstructure/exe:${PATH}"
ENV DATAPATH="/opt/RNAstructure/data_tables"

# LINEARFOLD (optional: the app falls back to windowed folding if it is missing)
RUN git clone --depth=1 \
    https://github.com/LinearFold/LinearFold.git \
    /opt/linearfold \
    && cd /opt/linearfold \
    && make \
    && chmod +x /opt/linearfold/linearfold \
    && echo "LinearFold installed" \
    || echo "WARNING: LinearFold build failed"

ENV PATH="/opt/linearfold:${PATH}"

# PYTHON APPLICATION
WORKDIR /app

COPY requirements.txt .
RUN pip3 install --no-cache-dir -r requirements.txt

COPY . .

# APPLICATION CONFIGURATION
ENV FLASK_DEBUG=False
ENV FLASK_PORT=8080
ENV PORT=8080
ENV DB_HOST=db
ENV DB_USER=root
ENV DB_NAME=cas13_grna_db
ENV MAX_SEQUENCE_LENGTH=20000
ENV GROQ_MODEL=openai/gpt-oss-20b

EXPOSE 8080

# STARTUP SCRIPT
COPY start.sh /start.sh

RUN sed -i 's/\r$//' /start.sh \
    && chmod +x /start.sh

CMD ["/start.sh"]
