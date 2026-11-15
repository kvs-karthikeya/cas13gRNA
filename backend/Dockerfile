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

# VIENNARNA
RUN RNAfold --version && echo "RNAfold OK"
RUN RNAplfold --version && echo "RNAplfold OK"

# RNASTRUCTURE 64-BIT TEXT INTERFACES
ARG RNASTRUCTURE_VERSION=6.6

RUN wget -q \
    "https://rna1.urmc.rochester.edu/Releases/current/RNAstructureLinuxTextInterfaces64bit.tgz" \
    -O /tmp/rnastructure.tgz \
    && tar -xzf /tmp/rnastructure.tgz -C /opt/ \
    && rm /tmp/rnastructure.tgz \
    && test -x /opt/RNAstructure/exe/Fold \
    && echo "RNAstructure Fold executable installed"

ENV PATH="/opt/RNAstructure/exe:${PATH}"
ENV DATAPATH="/opt/RNAstructure/data_tables"

# Verify RNAstructure
RUN command -v Fold && test -d "$DATAPATH"

# LINEARFOLD
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

