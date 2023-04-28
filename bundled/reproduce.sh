#!/bin/bash
# Reproduce results then text of the paper 
scriptDir="$(dirname "$0")"
cd "$(realpath $scriptDir)"

pip install -r ./binder/requirements.txt

# Regenerate computed results (figs) needed for compiling paper
./reproduce/computed.sh MIN # Replace with MAX to execute ApndxBalancedGrowthcNrmAndCov.ipynb
./reproduce/document.sh # Make latex document, figures, tables, etc
