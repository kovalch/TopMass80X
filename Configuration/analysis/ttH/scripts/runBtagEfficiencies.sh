#!/bin/sh


if [ $# -ge 1 ] ; then
    echo "No options allowed for this script. They are set within the script."
    echo "You used options: $@"
    exit 1
fi



source $(dirname `readlink -f $0`)/parallelTools.sh

systematics=(
    Nominal
    PU_UP PU_DOWN
    TRIG_UP TRIG_DOWN
    LEPT_UP LEPT_DOWN
    JER_UP JER_DOWN
    JES_UP JES_DOWN
    BTAG_UP BTAG_DOWN
    BTAG_LJET_UP BTAG_LJET_DOWN
    #KIN_UP KIN_DOWN
)

for syst in "${systematics[@]}" ; do
    for c in ee emu mumu ; do
        w
        $LA -f ttbarsignalplustau.root -p 0 -c $c -s $syst&
    done
done

wait

if [ "$isNAF" = 1 ]; then
    echo "Please check your jobs with qstat -u $USER | grep load_Analysis"
else
    echo "Processing all nominal samples finished!"
fi


