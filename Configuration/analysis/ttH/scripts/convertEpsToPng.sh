#! /bin/sh

# Run this script via (from TopAnalysis/Configuration/analysis/ttH) :
# ./scripts/convertEpsToPng.sh <path to folder with eps files>


InputFolder=$1

cd $InputFolder

for x in `ls *.eps`; do
y=`basename $x .eps`
convert $y.eps $y.png

done
