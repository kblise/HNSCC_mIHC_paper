# Single-Cell Spatial Proteomics Analyses of Head and Neck Squamous Cell Carcinoma Reveal Tumor Heterogeneity and Immune Architectures Associated with Clinical Outcome


This repository contains all code necessary to reproduce the results and figures of the paper "Single-Cell Spatial Proteomics Analyses of Head and Neck Squamous Cell Carcinoma Reveal Tumor Heterogeneity and Immune Architectures Associated with Clinical Outcome." All data, including the output of the multiplex immunohistochemistry computational imaging processing workflow for each tumor region and survival data for each patient are available on Zenodo: 10.5281/zenodo.4584441.


## Steps to run code

**0.** Install Conda: [Instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

**1.** Open terminal and navigate to desired directory: `cd Desired/Directory/Here/`

**2.** Download mIHC and clinical data files from Zenodo: `wget http://zenodo.org/record/4584441files/ZENODOFOLDERNAMEHERE.zip`

**3.** Unzip data folder: `unzip \*.zip`

**4.** Clone the repository: `git clone https://github.com/kblise/HNSCC_mIHC_paper.git`

**5.** Merge data folder with repo folder: `rsync -a ./ZENODOFOLDERHERE/ ./HNSCC_mIHC_paper/`

**6.** Navigate into HNSCC_mIHC_paper folder: `cd HNSCC_mIHC_paper/`

**7.** Create two new folders in HNSCC_mIHC_paper directory: `mkdir {dfCreated,figures}`

**8.** Create a new conda environment with the necessary packages: `conda env create -f hnscc_new_env.yml`

**9.** Activate the hnsccEnv environment: `source activate hnsccEnv` or `conda activate hnsccEnv`

**10.** From the command line, run python: `python`

**11.** Run hnsccMakeFigures.py program to create desired figure (1-4): 

`from hnsccMakeFigures import fig1`

`fig1()`



**Note: all csvs created to generate figures will be saved to the 'dfCreated' folder and all figures will be saved to the 'figures' folder.**
