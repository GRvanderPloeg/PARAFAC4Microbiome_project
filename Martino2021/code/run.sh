#!/usr/bin/env bash

: '
This is the master script for the capsule.
When you click "Reproducible Run", the code in this file will execute.
This script reproduces all of the figures from stored results. It
also includes some evaluations of results but does not include
the preprocessing from raw data. All of the resulting figures
and tables are output in the "results" directory.

The script called "run-complete.sh" follows the whole process from
raw data rto figure generation. However, it should be noted that running
the "run-complete.sh" script can take a very long time to complete
(up to several days). Therfore, it should likely be run locally 
or on a compute cluster.
'

timout_limit=1800
# PART 1: Introduction
simdir="1.0.0-simulation-benchmarking"
# simulation benchmarks
jupyter nbconvert --execute $simdir/2.2.0-Simulation-Benchmarking-Plotting.ipynb\
                  --ExecutePreprocessor.timeout=$timout_limit
# example IBD figure from Halfvarson et al.
jupyter nbconvert --execute $simdir/3.1.0-Benchmark-Plotting.ipynb\
                  --ExecutePreprocessor.timeout=$timout_limit

# PART 2: Case Studies
# two cases
case1="2.0.0-ECAM"
case2="2.1.0-DIABIMMUNE"
# both case studies are named the same
declare -a CaseNotebooks=("2.0-ordination-plotting.ipynb" "2.1-PERMAONVA.ipynb" "2.2-KNN.ipynb")
# Iterate the string array using for loop
for notebook in ${CaseNotebooks[@]}; do
   jupyter nbconvert --execute $case1/$notebook\
                  --ExecutePreprocessor.timeout=$timout_limit
   jupyter nbconvert --execute $case2/$notebook\
                  --ExecutePreprocessor.timeout=$timout_limit
done
# compare the feature ranks
compare="2.2.0-compare-feature-ranks"
jupyter nbconvert --execute $compare/1.0.0-rankings.ipynb\
                  --ExecutePreprocessor.timeout=$timout_limit

# PART 3: Meta-Analyses with log-ratios
# American Gut section
AG="3.0.0-AG-meta-analysis"
jupyter nbconvert --execute  $AG/1.0.0-AG-meta-analysis.ipynb\
                  --ExecutePreprocessor.timeout=$timout_limit

