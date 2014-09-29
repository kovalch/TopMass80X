#!/bin/bash

systematics=(
     Nominal
#    PU_UP PU_DOWN
#    TRIG_UP TRIG_DOWN
#    LEPT_UP LEPT_DOWN
#    JER_UP JER_DOWN
    JES_UP JES_DOWN
    BTAG_UP BTAG_DOWN
    BTAG_LJET_UP BTAG_LJET_DOWN
#    BTAGDISCR_BSTAT1_UP BTAGDISCR_BSTAT1_DOWN
#    BTAGDISCR_BSTAT2_UP BTAGDISCR_BSTAT2_DOWN
#    BTAGDISCR_LSTAT1_UP BTAGDISCR_LSTAT1_DOWN
#    BTAGDISCR_LSTAT2_UP BTAGDISCR_LSTAT2_DOWN
#    BTAGDISCR_CERR1_UP BTAGDISCR_CERR1_DOWN
#    BTAGDISCR_CERR2_UP BTAGDISCR_CERR2_DOWN
#    KIN_UP KIN_DOWN
#    TOP_PT_UP TOP_PT_DOWN
)

for systematic in ${systematics[@]} ; do
    echo "Producing input for systematic: "$systematic
    ./install/bin/Histo -d stacked -g dy ttbb -c combined -s ${systematic} 
done

echo "Production of all inputs finished"
echo
