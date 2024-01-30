#!/usr/bin/env bash

: '
This script runs all of the analyses for this method from
raw data to figure generation. All of the resulting figures
and tables are output in the "results" directory. However,
it should be noted that running the "run-complete.sh" script
can take a very long time to complete (up to several days).
Therfore, it should likely be run locally or on a compute cluster.
'

timout_limit=1800
# PART 1: Introduction
simdir="1.0.0-simulation-benchmarking"
# simulation benchmarks
declare -a SimNote=("1.0.0-dataset-preprocessing.ipynb" "2.1.0-QIIME2-Analysis-Per-Simulation.ipynb" "2.2.0-Simulation-Benchmarking-Plotting.ipynb")
for notebook in ${SimNote[@]}; do
   jupyter nbconvert --execute $simdir/$notebook\
                     --ExecutePreprocessor.timeout=$timout_limit
done
# example IBD figure from Halfvarson et al.
declare -a ExampleNote=("3.1.0-Benchmark-Plotting.ipynb" "3.0.0-Real-Data-Analysis.ipynb")
for notebook in ${ExampleNote[@]}; do
   jupyter nbconvert --execute $simdir/$notebook\
                     --ExecutePreprocessor.timeout=$timout_limit
done

# PART 2: Case Studies
# two cases
case1="2.0.0-ECAM"
case2="2.1.0-DIABIMMUNE"
# both case studies are named the same
declare -a CaseNotebooks=("1.0-preprocessing.ipynb" "1.1-QIIME2-analysis.ipynb" "2.0-ordination-plotting.ipynb" "2.1-PERMAONVA.ipynb" "2.2-KNN.ipynb")
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

