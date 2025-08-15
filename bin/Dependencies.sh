#!/bin/bash

# Make sure conda commands are available
source /home/$USER/miniconda3/etc/profile.d/conda.sh

# ---------------------------
# LTR_Finder 1.07 environment
# ---------------------------
conda create -n ltr_finder-1.07 python=3.10.6 -y
conda activate ltr_finder-1.07
conda install -c bioconda ltr_finder=1.07 -y
pip install pandas matplotlib biopython
conda deactivate

# ---------------------------
# LTR_Retriever 2.9.5 environment
# ---------------------------
conda create -n ltr_retriever-2.9.5 python=3.10.6 -y
conda activate ltr_retriever-2.9.5
conda install -c conda-forge LTR_retriever=2.9.5 -y
conda install pandas matplotlib -y
conda install -c bioconda mafft iqtree gffcompare -y
conda deactivate

# ---------------------------
# PANNZER2 environment
# ---------------------------
conda create -n pannzer2 python=3.7 -y
conda activate pannzer2
conda install -c conda-forge numpy scipy pandas fastcluster requests perl -y
conda deactivate
