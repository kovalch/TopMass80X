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
    powhegHerwig powheg.root mcatnlo
    Perugia11.root Perugia11NoCR
)


for systematic in "${systematics[@]}"; do
    printf "\n\n\n\033[1;1mStart running on systematic: ${systematic}\033[1;m\n\n\n"
    
    for channel in ee emu mumu; do
        for part in $(seq 0 4); do
            w
            $LA -f ttbarsignalplustau_${systematic} -p ${part} -c ${channel} $@ &
        done
    done
done


wait


if [ "$isNAF" = 1 ]; then
    echo "Please check your jobs with qstat -u $USER | grep load_Analysis"
else
    echo "Processing all theory variations finished!"
fi

