#! /bin/sh

InputFolder=$1

cd $InputFolder

for x in `ls *.eps`; do
y=`basename $x .eps`
convert $y.eps $y.png

done
