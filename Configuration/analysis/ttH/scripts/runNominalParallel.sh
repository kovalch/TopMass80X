#!/bin/sh


if [ $# -ge 1 ] && ( [[ "$@" == *"-c"* ]] || [[ "$@" == *"-s"* ]] || [[ "$@" == *"-f"* ]] || [[ "$@" == *"-p"* ]] ) ; then
    echo "Options '-s', '-c', '-f', '-p' are not allowed. They are set within the script."
    echo "You used options: $@"
    exit 1
fi



source $(dirname `readlink -f $0`)/parallelTools.sh


for channel in ee emu mumu; do
    w
    $LA -f ttbarsignalplustau.root -p 0 -c $channel $@ &
    $LA -f ttbarsignalplustau.root -p 101 -c $channel $@ &
    $LA -f ttbarsignalplustau.root -p 201 -c $channel $@ &
    w
    $LA -f ttbarsignalplustau.root -p 102 -c $channel $@ &
    $LA -f ttbarsignalplustau.root -p 202 -c $channel $@ &
    w
    $LA -f ttbarsignalplustau.root -p 103 -c $channel $@ &
    $LA -f ttbarsignalplustau.root -p 203 -c $channel $@ &
    $LA -f ttbarsignalplustau.root -p 4 -c $channel $@ &
    w
    $LA -f ttbarH125tobbbar -c $channel $@ &
    $LA -f ttbarH125incl -p 0 -c $channel $@ &
done

for channel in ee emu mumu; do
    w
    $LA -f dy50inf -p 0 -c $channel $@ &
    $LA -f dy50inf -p 1 -c $channel $@ &
    $LA -f dy50inf -p 2 -c $channel $@ &
done

for channel in ee emu mumu; do
    w
    $LA -f dy1050 -p 0 -c $channel $@ &
    $LA -f dy1050 -p 1 -c $channel $@ &
    $LA -f dy1050 -p 2 -c $channel $@ &
done

for channel in ee emu mumu; do
    w
    $LA -f ${c}_run2012A -c $channel $@ &
    $LA -f ${c}_run2012B -c $channel $@ &
    $LA -f ${c}_run2012C -c $channel $@ &
    $LA -f ${c}_run2012D -c $channel $@ &
done

for channel in ee emu mumu; do
    for pattern in qcd single ttbarbg.root wtol wwtoall wztoall zztoall ttbarW ttbarZ ttgamma www wwz zzz; do
        w
        $LA -f $pattern -c $channel $@ &
    done
done

wait

if [ "$isNAF" = 1 ]; then
    echo "Please check your jobs with qstat -u $USER | grep load_Analysis"
else
    echo "Processing all nominal samples finished!"
fi





