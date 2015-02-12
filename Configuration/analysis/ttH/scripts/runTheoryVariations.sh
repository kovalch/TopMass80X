#!/bin/sh

if [ $# -ge 1 ] && ( [[ "$@" == *"-c"* ]] || [[ "$@" == *"-s"* ]] || [[ "$@" == *"-f"* ]] || [[ "$@" == *"-p"* ]] ) ; then
    echo "Options '-s', '-c', '-f', '-p' are not allowed. They are set within the script."
    echo "You used options: $@"
    exit 1
fi

source $(dirname `readlink -f $0`)/parallelTools.sh

systematics=(
    massup massdown
    matchingup matchingdown
    scaleup scaledown
    powhegHerwig powheg mcatnlo
    Perugia11NoCR Perugia11
)


for systematic in "${systematics[@]}"; do
    printf "\n\n\n\033[1;1mStart running on systematic: ${systematic}\033[1;m\n\n\n"
    
    for channel in ee emu mumu; do
        w
        $LA -f ttbarsignalplustau_${systematic}.root -p 0 -c $channel $@ &
        $LA -f ttbarsignalplustau_${systematic}.root -p 101 -c $channel $@ &
        $LA -f ttbarsignalplustau_${systematic}.root -p 201 -c $channel $@ &
        w
        $LA -f ttbarsignalplustau_${systematic}.root -p 102 -c $channel $@ &
        $LA -f ttbarsignalplustau_${systematic}.root -p 202 -c $channel $@ &
        w
        $LA -f ttbarsignalplustau_${systematic}.root -p 103 -c $channel $@ &
        $LA -f ttbarsignalplustau_${systematic}.root -p 203 -c $channel $@ &
        $LA -f ttbarsignalplustau_${systematic}.root -p 4 -c $channel $@ &
    done
done


wait


if [ "$isNAF" = 1 ]; then
    echo "Please check your jobs with qstat -u $USER | grep load_Analysis"
else
    echo "Processing all theory variations finished!"
fi

