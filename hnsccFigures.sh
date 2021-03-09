#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code author: Katie E. Blise
Date: March, 2021

This .sh bash script will 
    -create a new conda environment
    -download clinical data and mIHC data from Zenodo DOI: 10.5281/zenodo.4584441
    -run hnsccMakeFigures.py file to generate Figures 1-4 from the following paper:
    'Single-cell spatial proteomics analyses of head and neck squamous cell carcinoma reveal tumor heterogeneity and immune architectures associated with clinical outcome'

"""

#Navigate into directory with github repo
cd HNSCC_mIHC_paper/

#Create two new folders in HNSCC_mIHC_paper directory
mkdir {dfCreated,figures}

#Create a new conda environment with necessary packages
conda env create -f hnscc_new_env.yml

#Acviate hnsccEnv conda environment
source activate hnsccEnv

#Install zenodo_get in hnsccEnv conda environmnet
pip install zenodo_get

#Download all data into HNSCC_mIHC_paper directory from Zenodo DOI: 10.5281/zenodo.4584441
zenodo_get -r 4584441

#Run python file hnsccMakeFigures.py to create Figures 1-4 
python
from hnsccMakeFigures import fig1, fig2, fig3, fig4
fig1()
fig2()
fig3()
fig4()

#Quit python
quit()

#Deactivate hnsccEnv conda environment
conda deactivate