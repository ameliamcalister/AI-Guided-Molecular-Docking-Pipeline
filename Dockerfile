# Ubuntu 24.04 ships both AutoDock Vina and Open Babel as apt packages, which
# the pipeline calls as command-line tools. They are the reason this runs as a
# container rather than on Streamlit Community Cloud.
#
# Reproducibility: the exact Vina version is recorded with every run by
# src/provenance, and a built image is immutable, so scores are reproducible
# for as long as a given image is deployed. Rebuilding may pick up a newer
# Vina from the archive; the build log below records which one was used.
FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      python3 python3-venv python3-pip openbabel autodock-vina ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# Recorded in the build log so a deployment can be traced to its tool versions.
RUN vina --version && obabel -V && python3 --version

RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

# Hugging Face Spaces runs the container as uid 1000, so everything the app
# writes to must be owned by that user. HOME matters because Streamlit writes
# caches there, and a read-only HOME fails at startup rather than at use.
#
# ubuntu:24.04 already ships a user "ubuntu" at uid 1000, so useradd would
# fail with "UID 1000 is not unique". Remove it first; the || true keeps this
# working on a base image that does not have it.
RUN (userdel -r ubuntu 2>/dev/null || true) \
 && useradd -m -u 1000 user
ENV HOME=/home/user
WORKDIR /app

# Run history and saved poses live here. On Spaces this is ephemeral: it is
# wiped whenever the Space restarts or rebuilds. Mount a volume to keep it:
#   docker run -v "$(pwd)/data:/app/data" ...
ENV RUNS_DB=/app/data/runs.db
RUN install -d -o user -g user /app/data

# --chown on each COPY rather than a recursive chown afterwards, which would
# duplicate every file into a new layer.
COPY --chown=user:user docking_pipeline.py docking_app_v2.py ./
COPY --chown=user:user src ./src
COPY --chown=user:user .streamlit ./.streamlit

USER user

EXPOSE 8501
CMD ["streamlit", "run", "docking_app_v2.py"]
