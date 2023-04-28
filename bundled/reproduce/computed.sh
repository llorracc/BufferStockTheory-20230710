#!/bin/bash

scriptDir="$(dirname "$0")" # Get the path to the directory this script is in
# scriptDir=/Volumes/Data/Papers/BufferStockTheory/BufferStockTheory-Latest/reproduce/
cd "$scriptDir/.."  # Move to its parent 

[[ ! -e ./binder/requirements.out ]] && echo ''  && echo 'Installing requirements' && echo '' && pip install -r ./binder/requirements.txt | tee binder/requirements.out 

echo '' ; echo 'Producing figures' ; echo ''

[[ ! -d ./Code/Python/src ]] && mkdir -p ./Code/Python/src

cp -r ./src/* ./Code/Python/src

cd "./Code/Python"
[[ ! -e BufferStockTheor*.py ]] && jupyter nbconvert --to script BufferStockTheor*.ipynb 
ipython BufferStockTheor*.py

rm BufferStockTheor*.py # Delete it to prevent jupytext conflicts

[[ -e latexdefs.tex ]] && rm -f latexdefs.tex # Delete junk file that might be created

cd ./Code/Python; ./test_Harmenbergs_method.sh

# Execute sims showing near-constant growth of mean c and cov(c,p), Ω_{M[c]} and Ω_{cov}
if [[ "$#" -gt 0 ]]; then
    if [[ "$1" != "MAX" ]]; then
	if [[ "$1" != "MIN" ]]; then
	    echo ''
	    echo "Only command line options are 'MIN' and 'MAX':"
	    echo ''
	    echo "./reproduce.sh MIN"
	    echo ''
	    echo 'Skips execution of the notebook Code/Python/ApndxBalancedGrowthcNrmAndCov.ipynb'
	    echo ''
	    echo "./reproduce.sh MAX"
	    echo ''
	    echo 'executes it.'
	    echo ''
	    echo 'That script requires large memory capacity, and takes'
	    echo 'many hours to run.  (You might want to do it overnight).'
	    echo ''
	else
	    if [[ "$1" == "MAX" ]]; then
		echo ipython ApndxBalancedGrowthcNrmAndCov.ipynb
	    fi
	fi
    fi
fi

# Delete the flag saying the requirements have been installed
# In case requirements have changed since last compilation
[[ -e binder/requirements.out ]] && rm -f binder/requirements.out
[[ -e src ]] && rm -Rf src

