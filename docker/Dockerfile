# SOLVIA custom ColabFold GPU image
# Base: NVIDIA CUDA Runtime for GPU support
FROM nvidia/cuda:12.2.0-runtime-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive

# System packages
RUN apt-get update && apt-get install -y \
    python3 python3-pip git wget curl build-essential \
    && rm -rf /var/lib/apt/lists/*

# Python upgrade
RUN python3 -m pip install --upgrade pip

# Core Python deps pinned for stability
RUN pip install \
    "colabfold==1.5.5" \
    "jax[cuda12_local]==0.4.31" \
    "jaxlib==0.4.31" \
    biopython==1.81 \
    numpy==1.24.3 \
    PyYAML==6.0.1

# Workdir
WORKDIR /work

# Default entrypoint runs colabfold_batch
ENTRYPOINT ["colabfold_batch"]
