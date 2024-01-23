#!/bin/sh

create_results=1    # run matlab results script (1) 


if [ $create_results -eq 1 ]
then

    cp ../src/params/create_param_files.m ../output/tables 
     
    pdflatex -shell-escape -output-directory=../output -halt-on-error ../output/results.tex | grep '^!.*' -A200 --color=always
    rm ../output/*.aux
    rm ../output/*.fdb_latexmk
    rm ../output/*.fls
    rm ../output/*.out
    rm ../output/*.log
fi


