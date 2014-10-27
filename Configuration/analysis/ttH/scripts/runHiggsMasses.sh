#!/bin/sh


if [ $# -ge 1 ] && ( [[ "$@" == *"-c"* ]] || [[ "$@" == *"-s"* ]] || [[ "$@" == *"-f"* ]] || [[ "$@" == *"-p"* ]] ) ; then
    echo "Options '-s', '-c', '-f', '-p' are not allowed. They are set within the script."
    echo "You used options: $@"
    exit 1
fi



source $(dirname `readlink -f $0`)/parallelTools.sh



for channel in ee emu mumu; do
    $LA -f ttbarH110tobbbar -c $channel $@ &
    $LA -f ttbarH110incl -c $channel -p 1 $@ &
    $LA -f ttbarH115tobbbar -c $channel $@ &
    w
    $LA -f ttbarH115incl -c $channel -p 1 $@ &
    $LA -f ttbarH120tobbbar -c $channel $@ &
    $LA -f ttbarH120incl -c $channel -p 1 $@ &
    w
    $LA -f ttbarH1225incl -c $channel -p 1 $@ &
    $LA -f ttbarH125tobbbar -c $channel $@ &
    $LA -f ttbarH125incl -c $channel -p 1 $@ &
    w
    $LA -f ttbarH1275incl -c $channel -p 1 $@ &
    $LA -f ttbarH130tobbbar -c $channel $@ &
    $LA -f ttbarH130incl -c $channel -p 1 $@ &
    w
    $LA -f ttbarH135tobbbar -c $channel $@ &
    $LA -f ttbarH135incl -c $channel -p 1 $@ &
    $LA -f ttbarH140incl -c $channel -p 1 $@ &
    w
done


wait

if [ "$isNAF" = 1 ]; then
    echo "Please check your jobs with qstat -u $USER | grep load_Analysis"
else
    echo "Processing all Higgs mass variation samples finished!"
fi





