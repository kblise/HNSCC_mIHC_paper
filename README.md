# Single-Cell Spatial Proteomics Analyses of Head and Neck Squamous Cell Carcinoma Reveal Tumor Heterogeneity and Immune Architectures Associated with Clinical Outcome


This repository contains all code necessary to reproduce the results and figures of the paper "Single-Cell Spatial Proteomics Analyses of Head and Neck Squamous Cell Carcinoma Reveal Tumor Heterogeneity and Immune Architectures Associated with Clinical Outcome." All data, including the output of the multiplex immunohistochemistry computational imaging processing workflow for each tumor region and survival data for each patient are available on Zenodo: [10.5281/zenodo.4584441](https://zenodo.org/record/4584441#.YEbMFZNKiYA).


## Steps to run code

**0.** Install Conda: [Instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

**1.** Open terminal and navigate to desired directory: `cd Desired/Directory/Here/`

**2.** Clone the repository: `git clone https://github.com/kblise/HNSCC_mIHC_paper.git`

**3.** Navigate into HNSCC_mIHC_paper folder: `cd HNSCC_mIHC_paper/`

**4.** Run the bash script to create Figures 1-4: `bash hnsccFigures.sh`


**Note: all csvs created to generate figures will be saved to the 'dfCreated' folder and all figures will be saved to the 'figures' folder.**


This program is intended for Python version 3.
