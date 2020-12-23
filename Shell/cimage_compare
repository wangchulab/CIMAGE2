#!/bin/bash

if [ $# -lt 3 ]; then
    echo Usage: $0 [by_protein] file1 column1 outname1 file2 column2 outname2 ...
    echo Align column1 in file1 with column2 in file2, and rename them as outname1 and outname2 in the output file
    exit -1
fi

if [ $1 == "by_protein" ]; then
    by_protein=$1
    shift
fi

cwd=$(pwd)


if [ "$by_protein" == "by_protein" ]; then
    R --vanilla --args $@ < $CIMAGE_PATH/R/compare_averaged_ratios_by_protein.R > compare_averaged_ratios_by_protein.Rout
else
    R --vanilla --args $@ < $CIMAGE_PATH/R/compare_averaged_ratios.R > compare_averaged_ratios.Rout
    outname=$(cat outname | sed 's/\.txt$//g')
    echo $cwd, $outname
    $CIMAGE_PATH/perl/textTableCompareToHtml.pl $outname $cwd
    add_preview $cwd/$outname.html
    rm -rf outname
fi
