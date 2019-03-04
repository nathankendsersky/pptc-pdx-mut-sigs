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
2. Code for Figure 4G

Details
=======

- deconstructsigs.sh
-- Run first to pull relevant files from figshare

- theme.R
-- Code for Figure 4G ggplot2 theme

- deconstructsigs-30-finalpub.R
-- Code to calculate cosine similarity value for each model and to plot Figure 4G
