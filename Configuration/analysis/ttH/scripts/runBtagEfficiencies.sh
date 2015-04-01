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
    KIN_UP KIN_DOWN
    TOP_PT_UP TOP_PT_DOWN
)

for systematic in "${systematics[@]}" ; do
    printf "\n\n\n\033[1;1mStart running on systematic: ${systematic}\033[1;m\n\n\n"
    
    for channel in ee emu mumu ; do
        w
        $LA -f ttbarsignalplustau.root -p 0 -c $channel -s $systematic &
    done
done

wait

if [ "$isNAF" = 1 ]; then
    echo "Please check your jobs with qstat -u $USER | grep load_Analy"
else
    echo "Processing all samples for b-tag efficiencies finished!"
fi


