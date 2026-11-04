# Ubuntu 24.04 ships both AutoDock Vina (1.2.5) and OpenBabel as apt packages,
# which the pipeline calls as command-line tools.
FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive PYTHONUNBUFFERED=1

RUN apt-get update \
 && apt-get install -y --no-install-recommends python3 python3-venv python3-pip openbabel autodock-vina \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /app
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY docking_pipeline.py docking_app_v2.py ./
COPY src ./src

# Run history and saved poses live here. Mount a folder to keep them: -v "$(pwd)/data:/app/data"
ENV RUNS_DB=/app/data/runs.db
RUN mkdir -p /app/data

EXPOSE 8501
CMD ["streamlit", "run", "docking_app_v2.py", "--server.address=0.0.0.0", "--server.port=8501"]
