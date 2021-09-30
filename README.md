# Single-Cell Spatial Proteomics Analyses of Head and Neck Squamous Cell Carcinoma Reveal Tumor Heterogeneity and Immune Architectures Associated with Clinical Outcome


This repository contains all code necessary to reproduce the results and figures of the paper "Single-Cell Spatial Proteomics Analyses of Head and Neck Squamous Cell Carcinoma Reveal Tumor Heterogeneity and Immune Architectures Associated with Clinical Outcome." All data, including the output of the multiplex immunohistochemistry computational imaging processing workflow for each tumor region, and clinical data for each patient are available on Zenodo: [DOI: 10.5281/zenodo.5540356](https://doi.org/10.5281/zenodo.5540356).

# Steps to Create Figures

## 1) Clone repository

**a.** Open Terminal and navigate to desired directory: `cd Desired/Directory/Here/`

**b.** Clone repo: `git clone https://github.com/kblise/HNSCC_mIHC_paper.git`

## 2) Create new conda environment

**a.** Install Conda if not already installed: [Instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

**b.** Navigate into HNSCC_mIHC_paper directory: `cd HNSCC_mIHC_paper/`

**b.** Create new conda environment with necessary packages: `conda env create -f hnscc_new_env.yml`

**c.** Activate hnsccEnv conda environment: `source activate hnsccEnv` or `conda activate hnsccEnv`

## 3) Create new folders and download data from Zenodo

**a.** Create three new folders in HNSCC_mIHC_paper directory: `mkdir {data,dfCreated,figures}`

**b.** Navigate into data directory: `cd data/`

**c.** Install zenodo_get in hnsccEnv conda environment: `pip install zenodo_get`

**d.** Download all data into HNSCC_mIHC_paper/data/ directory from Zenodo: `zenodo_get -r 5540356`

## 4) Create Figures 1-4, Table 3, and Supplemental Figures 1-3

**a.** Navigate back to HNSCC_mIHC_paper directory: `cd ..`

**b.** Run python file hnsccMakeFigures.py to create Figures 1-4, Table 3, and Supplemental Figures 1-3: `python hnsccMakeFigures.py`



**Note: all csvs created to generate figures will be saved to the 'dfCreated' folder and all figures and tables will be saved to the 'figures' folder.**


This program is intended for Python version 3.
