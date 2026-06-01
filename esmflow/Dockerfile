FROM nvidia/cuda:12.4.1-cudnn-devel-ubuntu22.04
LABEL description="Optimized environment for running AlphaFlow (PR 58) with CUDA 12."

# 1. Install critical system utilities
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update -y && apt-get install -y --no-install-recommends \
        git \
        tmux \
        wget \
        ca-certificates \
        ninja-build \
    && rm -rf /var/lib/apt/lists/*

# 2. Define global environment working directories
ENV WORKSPACE_DIR="/opt/workspace"
ENV CONDA_DIR="/opt/mambaforge"
ENV PATH="$CONDA_DIR/bin:$PATH"
WORKDIR $WORKSPACE_DIR

# 3. Install Miniforge 
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O miniforge.sh \
    && chmod +x miniforge.sh \
    && ./miniforge.sh -b -p $CONDA_DIR \
    && rm miniforge.sh

# 4. Create Conda environment with Python 3.9
RUN mamba create -n alphaflow python=3.9 -y

# Use conda environment shell for subsequent setup steps
SHELL ["conda", "run", "-n", "alphaflow", "/bin/bash", "-c"]

### Start of AlphaFlow specific part ###

# 5. Install PyTorch 2.5.1 (defaults to CUDA 12.4 compatibility)
RUN pip install --no-cache-dir torch==2.5.1 torchvision torchaudio

# 6. Install base utilities and explicit pinned dependencies
RUN pip install --no-cache-dir tqdm numpy==1.21.2 pandas==1.5.3 scipy==1.7.1

# 7. Install structural dependencies via pip to keep NumPy bindings aligned
RUN pip install --no-cache-dir biopython dm-tree==0.1.6 modelcif==0.7 ml-collections==0.1.0 absl-py einops
RUN mamba install -c conda-forge mdtraj==1.9.9 -y

# 8. Install deep learning wrappers (Removed redundant torch installation)
RUN pip install --no-cache-dir "pytorch_lightning==2.0.4" "fair-esm" "wandb"

# 9. Install OpenFold from the pl_upgrades branch for PyTorch 2.x/CUDA 12 support
ENV FORCE_CUDA=1
RUN pip install --no-cache-dir 'openfold @ git+https://github.com/aqlaboratory/openfold.git@pl_upgrades'

# 10. Clone the AlphaFlow pipeline code repository and checkout PR 58
RUN git clone https://github.com/bjing2016/alphaflow.git $WORKSPACE_DIR/alphaflow \
    && cd $WORKSPACE_DIR/alphaflow \
    && git fetch origin pull/58/head:pr58 \
    && git checkout pr58

# 11. Download ESMFlow model weights during image build
RUN mkdir -p $WORKSPACE_DIR/alphaflow/params \
    && wget --progress=dot:giga https://huggingface.co/bjing-mit/alphaflow/resolve/main/params/esmflow_pdb_base_202402.pt \
       -O $WORKSPACE_DIR/alphaflow/params/esmflow_pdb_base_202402.pt

### End of AlphaFlow specific part ###

# 12. Default initialization setups
SHELL ["/bin/bash", "-c"]
RUN conda init bash \
    && echo "conda activate alphaflow" >> ~/.bashrc

ENTRYPOINT ["/bin/bash"]