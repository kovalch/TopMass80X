#!/bin/sh


if [ $# -ge 1 ] && ( [[ "$@" == *"-c"* ]] || [[ "$@" == *"-s"* ]] || [[ "$@" == *"-f"* ]] || [[ "$@" == *"-p"* ]] ) ; then
    echo "Options '-s', '-c', '-f', '-p' are not allowed. They are set within the script."
    echo "You used options: $@"
    exit 1
fi



source $(dirname `readlink -f $0`)/parallelTools.sh


for c in ee emu mumu; do
    w
    $LA -f dy -p 0 -c $c $@ &
    $LA -f dy -p 1 -c $c $@ &
    $LA -f dy -p 2 -c $c $@ &
    $LA -f ttbarsignalplustau.root -p 0 -c $c $@ &
    $LA -f ttbarsignalplustau.root -p 1 -c $c $@ &
    $LA -f ttbarsignalplustau.root -p 2 -c $c $@ &
    $LA -f ttbarH125tobbbar -c $c $@ &
    $LA -f ttbarH125incl -p 0 -c $c $@ &
done

for c in ee emu mumu; do
    w
    $LA -f ${c}_run2012A -c $c $@ &
    $LA -f ${c}_run2012B -c $c $@ &
    $LA -f ${c}_run2012C -c $c $@ &
    $LA -f ${c}_run2012D -c $c $@ &
done

for i in qcd single ttbarbg.root wtol wwtoall wztoall zztoall ttbarW ttbarZ; do
    w
    $LA -f $i -c ee $@ &
    $LA -f $i -c emu $@ &
    $LA -f $i -c mumu $@ &
done

wait

if [ "$isNAF" = 1 ]; then
    echo "Please check your jobs with qstat -u $USER | grep load_Analysis"
else
    echo "Processing all nominal samples finished!"
fi





