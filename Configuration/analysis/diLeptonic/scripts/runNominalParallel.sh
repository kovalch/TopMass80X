#!/bin/sh

#
# $LA -f dy -d 11 -c ee
# $LA -f dy -d 13 -c mumu
# $LA -f ee_run2012B -c ee
# $LA -f ee_run2012C -c ee
# $LA -f mumu_run2012B -c mumu
# $LA -f mumu_run2012C -c mumu
#

source $(dirname `readlink -f $0`)/parallelTools.sh

for c in ee emu mumu; do
    w
    $LA -f ttbarsignalplustau.root -c $c &
    $LA -f dy -d 11 -c $c &
    $LA -f dy -d 13 -c $c &
    $LA -f dy -d 15 -c $c &
done

for c in ee emu mumu; do
    w
    $LA -f ${c}_run2012A -c $c &
    $LA -f ${c}_run2012B -c $c &
    $LA -f ${c}_run2012C -c $c &
    $LA -f ${c}_run2012D -c $c &
done 

for i in zz qcd single ttbarbg.root wtol ww wz ttbarW ttgjets ttbarZ; do
    w
    if [ "$isNAF" = 1 ]; then
        $LA -f $i
    else
        $LA -f $i -c ee&
        $LA -f $i -c emu&
        $LA -f $i -c mumu&
    fi
done

wait

if [ "$isNAF" = 1 ]; then
    echo "Please check your jobs with qstat -u $USER | grep load_Analy"
else
    echo "Processing all nominal samples finished!"
fi


