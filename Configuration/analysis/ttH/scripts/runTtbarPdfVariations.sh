#!/bin/sh

if [ $# -ge 1 ] && ( [[ "$@" == *"-c"* ]] || [[ "$@" == *"-s"* ]] || [[ "$@" == *"-f"* ]] || [[ "$@" == *"-p"* ]] ) ; then
    echo "Options '-s', '-c', '-f', '-p' are not allowed. They are set within the script."
    echo "You used options: $@"
    exit 1
fi

source $(dirname `readlink -f $0`)/parallelTools.sh


for number in `seq 0 52`; do
    printf "\n\n\n\033[1;1mStart running PDF variation number: ${number}\033[1;m\n\n\n"
    
    for channel in ee emu mumu; do
        w
        $LA -f ttbarsignalplustau.root -p 0 -c $channel -s PDF --sid $number $@ &
        $LA -f ttbarsignalplustau.root -p 101 -c $channel -s PDF --sid $number $@ &
        $LA -f ttbarsignalplustau.root -p 201 -c $channel -s PDF --sid $number $@ &
        w
        $LA -f ttbarsignalplustau.root -p 102 -c $channel -s PDF --sid $number $@ &
        $LA -f ttbarsignalplustau.root -p 202 -c $channel -s PDF --sid $number $@ &
        w
        $LA -f ttbarsignalplustau.root -p 103 -c $channel -s PDF --sid $number $@ &
        $LA -f ttbarsignalplustau.root -p 203 -c $channel -s PDF --sid $number $@ &
        $LA -f ttbarsignalplustau.root -p 4 -c $channel -s PDF --sid $number $@ &
    done
done


wait


if [ "$isNAF" = 1 ]; then
    echo "Please check your jobs with qstat -u $USER | grep load_Analy"
else
    echo "Processing all ttbar PDF variations finished!"
fi

