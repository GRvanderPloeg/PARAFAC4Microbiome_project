# Benchmarking and Analysis for CTF 

CTF is a proposed method for running tensor factorization on sparse compositional omics data sets. [Gemelli](https://github.com/cameronmartino/gemelli/tree/master/gemelli). (the CTF codebase) performs unsupervised dimensionality reduction of spatiotemporal microbiome data. The output of gemelli helps to resolve spatiotemporal subject variation and the biological features that separate them. 

Here we provide all of the code used to perform benchmarking and analyze the case studies performed for the CTF publication.

## Installation 

Gemelli is most easily used inside of a [QIIME2](https://qiime2.org/) environment. The directions for creating a QIIME2 environment can be found [here](https://docs.qiime2.org/2019.10/install/native/#install-qiime-2-within-a-conda-environment). All of the versioning is stored in the dockerfile in this capsule.

# Run files

There are two options for running the analyses. The less compute intensive option `run.sh` is to only re-generate the figures from the processing outputs. The other option `run-complete.sh` is to run everything from data acquired from their respective Qiita processing repositories (outlined below). 

# Data

All of the processed data sets used in this project are open source. Below is an outline of the data sets used with links to acquire them from the source. A processing version is also stored in the data directory of this capsule. 

* data driven simulation benchmarking
    * `Halfvarson-IBD-Qiita-1629`: Qiita study ID [1629](https://qiita.ucsd.edu/study/description/1629)
* case studies in infant development
    * `ECAM-Qiita-10249`: Qiita study ID [10249](https://qiita.ucsd.edu/study/description/10249)
    * `DIABIMMUNE-Qiita-11884`: Qiita study ID [11884](https://qiita.ucsd.edu/study/description/11884)
* meta analysis using shared ranks from case studies
    * American Gut [ftp.microbio.me](ftp://ftp.microbio.me/AmericanGut/manuscript-package/2500/)

# Code

## (1) Data-Driven Simulation Benchmarking

All of the simulation benchmarking code is stored in `1.0.0-simulation-benchmarking`. The jupyter notebooks it contains are numerically ordered by the order they need to be run. These notebooks have the helper functions in `simulation-code` for performing the data-driven simulation benchmarking based on the simulation used on [Äijö et al](https://www.ncbi.nlm.nih.gov/pubmed/28968799) and run on the IBD data set acquired from Qiita study ID [1629](https://qiita.ucsd.edu/study/description/1629).

* 1.0.0-simulation-benchmarking
    * 1.0.0-dataset-preprocessing.ipynb
    * 2.1.0-QIIME2-Analysis-Per-Simulation.ipynb
    * 2.2.0-Simulation-Benchmarking-Plotting.ipynb
    * 3.0.0-Real-Data-Analysis.ipynb
    * 3.1.0-Benchmark-Plotting.ipynb

## (2) Infant Development Case Studies: ECAM & DIABIMMUNE

In order to benchmark CTF we used to real data case studies and compared the results to the commonly used existing methodology for dimensionality reduction. All of the data used for the case studies is publicly available on Qiita. For the first case study, the ECAM data set published by [Bokulich et al](https://www.ncbi.nlm.nih.gov/pubmed/27306664) the data is available with the Qiita ID [10249](https://qiita.ucsd.edu/study/description/10249). For the second case study, the DIABIMMUNE data set published by [Yassour et al](https://stm.sciencemag.org/content/8/343/343ra81) the data is available with the Qiita ID [11884](https://qiita.ucsd.edu/study/description/11884). Both of the data sets analyses are structured as data preprocessing (1.0), QIIME2 processing for core tools as well as CTF (1.1), ordination plotting and log ratios (2.0), PERMANOVA (2.1), and finally KNN-classification APR (2.2). This follows the analysis from raw data tables acquired from Qiita to all benchmarking figure generation.

* 2.0.0-ECAM
    * 1.0-preprocessing.ipynb
    * 1.1-QIIME2-analysis.ipynb
    * 2.0-ordination-plotting.ipynb
    * 2.1-PERMAONVA.ipynb
    * 2.2-KNN.ipynb
* 2.1.0-DIABIMMUNE
    * 1.0-preprocessing.ipynb
    * 1.1-QIIME2-analysis.ipynb
    * 2.0-ordination-plotting.ipynb
    * 2.1-PERMAONVA.ipynb
    * 2.2-KNN.ipynb

Next we compare the feature ranking between the two case study data sets to evaluate the reproducibility between independent data sets.

* 2.2.0-compare-feature-ranks
    * 1.0.0-rankings.ipynb

## (3) Meta-analysis from feature rankings

First, we performed meta-analysis on the American Gut data set using the shared rankings from the ECAM and DIABIMMUNE results. The data acquisition results for this analysis is provided in a single notebook located in:

* 3.0.0-AG-meta-analysis/1.0.0-AG-meta-analysis.ipynb
