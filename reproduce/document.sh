#!/bin/bash

echo '' ; echo 'Reproduce text of paper:' ; echo ''

DIR=$(dirname $(realpath $0))

# Make sure tlmgr (texlive manager) is installed and initialized
[[ "$(which tlmgr)" == "" ]] && echo 'tlmgr is not available; install texlive and rerun'
[[ "$(which ~/.tlpkg)" == "" ]] && tlmgr init-usertree
# DIR=/Volumes/Data/Papers/BufferStockTheory/BufferStockTheory-Latest/reproduce/
texname=BufferStockTheory
output_directory='LaTeX'

cd "$DIR/.."

# Make figures that get made by executing a latex file
# (they should have a filename ending in Make.tex)
cd Figures
# For this paper, only the tikz figures need to be made by pdflatex - others are made by python
for fName_tikzMake in *Make.tex; do # names of all files ending in Make.tex
    echo "Processing figure $fName_tikzMake"
    fName=${fName_tikzMake%_tikzMake.tex} # Remove the '_tikzMake.tex' part of the filename
#    dep="pwd ; texliveonfly $fName_tikzMake"
#    echo dep="$dep"
#    eval "$dep"
    cmd="pdflatex -halt-on-error --output-format pdf -output-directory=../$output_directory $fName_tikzMake"
    echo "$cmd"
    eval "$cmd"
    echo    mv -f                                                             "../$output_directory/$fName.pdf" "$fName.pdf"
    mv -f                                                             "../$output_directory/$fName"_tikzMake".pdf" "$fName.pdf"    
done
cd ..

# Compile LaTeX files in root directory
for file in "$texname" "$texname"-NoAppendix "$texname"-Slides; do
    echo '' ; echo "Compiling $file" ; echo ''
#    dep="pwd ; texliveonfly $file"
#    echo dep="$dep"
#    eval "$dep"
    cmd="pdflatex -halt-on-error -output-directory=$output_directory $file"
    eval "$cmd"
    eval "$cmd > /dev/null" # Hide second output to reduce clutter
    bibtex $output_directory/"$file"
    eval "$cmd" # Hide third output to reduce clutter
    eval "$cmd > /dev/null" 
    echo '' ; echo "Compiled $file" ; echo ''
done

# Compile All-Figures and All-Tables
for type in Figures Tables; do
    # dep="texliveonfly $type/All-$type"
    # echo "pwd ; $dep"
    # eval "$dep"
    cmd="pdflatex -halt-on-error -output-directory=$output_directory $type/All-$type"
    echo "$cmd" ; eval "$cmd"
    # If there is a .bib file, make the references
    [[ -e "../$output_directory/$type/All-$type.aux" ]] && bibtex "$type/All-$type.bib" && eval "$cmd" && eval "$cmd" 
    mv -f "$output_directory/All-$type.pdf" "$type"  # Move from the LaTeX output directory to the destination
done

# All the appendices can be compiled as standalone documents (they are "subfiles")
# Make a list of all the appendices, put the list in the file /tmp/appendices
find ./Appendices -name '*.tex' ! -name '*econtexRoot*' ! -name '*econtexPath*' -maxdepth 1 -exec echo {} \; > /tmp/appendices

# For each appendix process it by pdflatex
# If it contains a standalone bibliography, process that
# Then rerun pdflatex to complete the processing and move the resulting pdf file

while read appendixName; do
    filename=$(basename ${appendixName%.*}) # Strip the path and the ".tex"
#    dep="texliveonfly $filename"
#    echo dep="$dep"
#    eval "$dep"
    cmd="pdflatex -halt-on-error                 --output-directory=$output_directory $appendixName"
    echo "$cmd"
    eval "$cmd"
    if grep -q 'bibliography{' "$appendixName"; then # it has a bibliography
	bibtex $output_directory/$filename 
	eval "$cmd" 
    fi
    eval "$cmd"
    cmd="mv $output_directory/$filename.pdf Appendices"
    echo "$cmd"
    eval "$cmd"
done < /tmp/appendices

[[ -e "$texname".pdf ]] && rm -f "$texname".pdf

echo '' 

if [[ -e "$output_directory/$texname.pdf" ]]; then
    echo "Paper has been compiled to $output_directory/$texname.pdf"
    echo "and copied to ./$texname.pdf"
    cp "$output_directory/$texname.pdf" "./$texname.pdf"
else
    echo "Something went wrong and the paper is not in $output_directory/$texname.pdf"
fi

echo ''

