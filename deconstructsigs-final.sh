#!/bin/bash

####################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create Directories
## change destination of folder in in lines 7 & 8, if desired
mkdir ~/Desktop/mutsig-demo/
cd ~/Desktop/mutsig-demo/
mkdir ./signatures
mkdir ./signatures/bysample
mkdir ./figures
mkdir ./figures/signatures

####################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Download MAF file
# https://figshare.com/projects/Genomic_landscape_of_childhood_cancer_patient-derived_xenograft_models/38147
# download requires wget (homebrew)
# $ brew install wget

# ##################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MAF file
wget --output-document='2019-02-14-allpdx-clean-maf-240.rda' https://ndownloader.figshare.com/files/14414198

# ##################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Clinical file
wget --output-document='pptc-pdx-clinical-web.txt' https://ndownloader.figshare.com/files/16603061

# ##################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Muts/Mb file
wget --output-document='2019-03-14-mutations-per-model.txt' https://ndownloader.figshare.com/files/15597293

# ##################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ggplot theme
wget --output-document='theme.R' https://raw.githubusercontent.com/marislab/pptc-pdx-mut-sigs/master/theme.R

# if wget is NOT installed...
# 1) install homebrew (https://brew.sh)
# 2) install wget (https://mikebeach.org/2012/09/29/how-to-install-wget-in-mac-os-x/)