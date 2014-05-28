#!/bin/sh

#######################################################
###  This scripts creates the BTagging efficiency 
###  histograms for each systematic variation
#######################################################

source $(dirname `readlink -f $0`)/parallelTools.sh

for syst in Nominal \
            JER_UP JER_DOWN \
            JES_UP JES_DOWN \
            PU_UP PU_DOWN \
            TRIG_UP TRIG_DOWN \
            LEPT_UP LEPT_DOWN \
            KIN_UP KIN_DOWN \
            BTAG_UP BTAG_DOWN BTAG_LJET_UP BTAG_LJET_DOWN \
            BTAG_PT_UP BTAG_PT_DOWN BTAG_ETA_UP BTAG_ETA_DOWN \
            BTAG_LJET_PT_UP BTAG_LJET_PT_DOWN BTAG_LJET_ETA_UP BTAG_LJET_ETA_DOWN; do
    for c in ee emu mumu; do
        w
        $LA -f ttbarsignalplustau.root -c $c -s $syst&
    done
done

#######################################################
###  Signal variations
###  PDF systematic will use the Nominal btag effciencies
#######################################################
for c in ee emu mumu; do
    w
    $LA -f tau_massup.root -c $c &
    $LA -f tau_massdown.root -c $c &
    $LA -f tau_matchingup.root -c $c &
    $LA -f tau_matchingdown.root -c $c &
    $LA -f tau_scaleup.root -c $c &
    $LA -f tau_scaledown.root -c $c &
    $LA -f tau_mcatnlo.root -c $c &
    $LA -f tau_powheg.root -c $c &
    $LA -f tau_powhegHerwig.root -c $c &
done



wait

if [ "$isNAF" = 1 ]; then
    echo "Please check your jobs with qstat -u $USER | grep load_Analysis"
else
    echo "Processing all nominal samples finished!"
fi


