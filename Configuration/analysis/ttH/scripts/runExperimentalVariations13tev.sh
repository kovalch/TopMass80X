#!/bin/sh


if [ $# -ge 1 ] && ( [[ "$@" == *"-c"* ]] || [[ "$@" == *"-s"* ]] || [[ "$@" == *"-f"* ]] || [[ "$@" == *"-p"* ]] ) ; then
    echo "Options '-s', '-c', '-f', '-p' are not allowed. They are set within the script."
    echo "You used options: $@"
    exit 1
fi



source $(dirname `readlink -f $0`)/parallelTools.sh


systematics=(
#     PU_UP PU_DOWN
#     TRIG_UP TRIG_DOWN
#     LEPT_UP LEPT_DOWN
#     JER_UP JER_DOWN
    JES_UP JES_DOWN
    BTAGDISCR_BPURITY_UP BTAGDISCR_BPURITY_DOWN
    BTAGDISCR_LPURITY_UP BTAGDISCR_LPURITY_DOWN
    BTAGDISCR_BSTAT1_UP BTAGDISCR_BSTAT1_DOWN
    BTAGDISCR_BSTAT2_UP BTAGDISCR_BSTAT2_DOWN
    BTAGDISCR_LSTAT1_UP BTAGDISCR_LSTAT1_DOWN
    BTAGDISCR_LSTAT2_UP BTAGDISCR_LSTAT2_DOWN
    BTAGDISCR_CERR1_UP BTAGDISCR_CERR1_DOWN
    BTAGDISCR_CERR2_UP BTAGDISCR_CERR2_DOWN
#     KIN_UP KIN_DOWN
#     TOP_PT_UP TOP_PT_DOWN
)


for systematic in "${systematics[@]}" ; do
    printf "\n\n\n\033[1;1mStart running on systematic: ${systematic}\033[1;m\n\n\n"
    
    for channel in ee emu mumu; do
        w
        $LA -f ttbarsignalplustau.root -p 0 -c $channel -s $systematic $@ &
        $LA -f ttbarsignalplustau.root -p 101 -c $channel -s $systematic $@ &
        $LA -f ttbarsignalplustau.root -p 201 -c $channel -s $systematic $@ &
        w
        $LA -f ttbarsignalplustau.root -p 102 -c $channel -s $systematic $@ &
        $LA -f ttbarsignalplustau.root -p 202 -c $channel -s $systematic $@ &
        w
        $LA -f ttbarsignalplustau.root -p 103 -c $channel -s $systematic $@ &
        $LA -f ttbarsignalplustau.root -p 203 -c $channel -s $systematic $@ &
        $LA -f ttbarsignalplustau.root -p 4 -c $channel -s $systematic $@ &
        w
        $LA -f ttbarH125tobbbar -c $channel -s $systematic $@ &
        $LA -f ttbarH125incl -p 0 -c $channel -s $systematic $@ &
        $LA -f ttbarH125incl -p 1 -c $channel -s $systematic $@ &
    done
    
    for channel in ee emu mumu; do
        w
        $LA -f dy50inf -p 0 -c $channel -s $systematic $@ &
        $LA -f dy50inf -p 1 -c $channel -s $systematic $@ &
        $LA -f dy50inf -p 2 -c $channel -s $systematic $@ &
    done
    
    for channel in ee emu mumu; do
        w
        $LA -f dy1050 -p 0 -c $channel -s $systematic $@ &
        $LA -f dy1050 -p 1 -c $channel -s $systematic $@ &
        $LA -f dy1050 -p 2 -c $channel -s $systematic $@ &
    done
    
    for channel in ee emu mumu; do
        w
        $LA -f ttbarbg.root -p 0 -c $channel -s $systematic $@ &
        $LA -f ttbarbg.root -p 1 -c $channel -s $systematic $@ &
        $LA -f ttbarbg.root -p 2 -c $channel -s $systematic $@ &
        w
        $LA -f ttbarbg.root -p 3 -c $channel -s $systematic $@ &
        $LA -f ttbarbg.root -p 4 -c $channel -s $systematic $@ &
    done

    for channel in ee emu mumu; do
        for pattern in qcd single wtol wwtoall wztoall zztoall ttbarW ttbarZ ttgamma www wwz zzz; do
            w
            $LA -f $pattern -c $channel -s $systematic $@ &
        done
    done
done


wait

if [ "$isNAF" = 1 ]; then
    echo "Please check your jobs with qstat -u $USER | grep load_Analy"
else
    echo "Processing all variations of nominal samples finished!"
fi





