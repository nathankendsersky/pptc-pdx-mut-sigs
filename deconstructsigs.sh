#!/bin/bash

####################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Download MAF file
# https://figshare.com/projects/Genomic_landscape_of_childhood_cancer_patient-derived_xenograft_models/38147


# ##################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MAF file
wget --output-document='2019-02-14-allpdx-clean-maf-240.rda' https://ndownloader.figshare.com/files/14414198

# ##################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Clinical file
wget --output-document='pptc-pdx-clinical-web.txt' https://ndownloader.figshare.com/files/14508536