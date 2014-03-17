#!/bin/sh

source $(dirname `readlink -f $0`)/parallelTools.sh

for sys in JES_UP JES_DOWN JER_UP JER_DOWN \
           PU_UP PU_DOWN \
           TRIG_UP TRIG_DOWN \
           LEPT_UP LEPT_DOWN \
           KIN_UP KIN_DOWN \
           BTAG_UP BTAG_DOWN BTAG_LJET_UP BTAG_LJET_DOWN \
           BTAG_PT_UP BTAG_PT_DOWN BTAG_ETA_UP BTAG_ETA_DOWN \
           BTAG_LJET_PT_UP BTAG_LJET_PT_DOWN BTAG_LJET_ETA_UP BTAG_LJET_ETA_DOWN; do # \
#            BTAG_BEFF_UP BTAG_BEFF_DOWN BTAG_CEFF_UP BTAG_CEFF_DOWN BTAG_LEFF_UP BTAG_LEFF_DOWN; do

    for c in ee emu mumu; do
        w
        $LA -f ttbarsignalplustau.root -c $c -s $sys&
        $LA -f dy -d 11 -c $c -s $sys&
        $LA -f dy -d 13 -c $c -s $sys&
        $LA -f dy -d 15 -c $c -s $sys&
    done

    for i in qcd single ttbarbg.root wtol ww wz zz ttbarW ttgjets ttbarZ; do
        w
        $LA -f $i -s $sys&
    done

done
wait

