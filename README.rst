.. |date| date::

*******************************
PPTC PDX Mutational Signatures
*******************************

:authors: Jo Lynne Rokita, Nathan Kendsersky
:contact: Jo Lynne Rokita (rokita@email.chop.edu)
:organization: CHOP
:status: In-process
:date: |date|

.. meta::
   :keywords: pdx, mouse, WES, COSMIC, mutational signatures, 2019
   :description: pdx WES somatic mutational signature analysis

Introduction
============

This repo contains code for:

1. Creation of mutational signature weight matrices
2. Code for Figure S3A

Details
=======
- Download the two following scripts.
  * deconstructsigs-final.sh -- RUN 1st -- pulls relevant files from figshare and create appropriate directories
  * deconstructsigs-final.R -- RUN 2nd -- calculates cosine similarity values and plot signatures by model/histology
- The following files will be loaded from the previous two scripts - no need to download.
  * theme.R -- contains theme for ggplot
  * leukemia-brain-order.txt -- contains the pdx order for ggplot (leukemias and brain tumors)
  * solid-order.txt -- contains the pdx order for ggplot (solid tumors)
  * signature-category-color.txt -- contains signatures, categories, and colors (for ggplot)
