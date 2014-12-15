#!/bin/bash

systematics=(
    Nominal
    JES_UP JES_DOWN
    BTAGDISCR_PURITY_UP BTAGDISCR_PURITY_DOWN
    BTAGDISCR_BSTAT1_UP BTAGDISCR_BSTAT1_DOWN
    BTAGDISCR_BSTAT2_UP BTAGDISCR_BSTAT2_DOWN
    BTAGDISCR_LSTAT1_UP BTAGDISCR_LSTAT1_DOWN
    BTAGDISCR_LSTAT2_UP BTAGDISCR_LSTAT2_DOWN
    XSEC_TT2B_UP XSEC_TT2B_DOWN
    XSEC_TTCC_UP XSEC_TTCC_DOWN
#    PU_UP PU_DOWN
#    TRIG_UP TRIG_DOWN
#    LEPT_UP LEPT_DOWN
#    JER_UP JER_DOWN
#    KIN_UP KIN_DOWN
#    TOP_PT_UP TOP_PT_DOWN
)

for systematic in ${systematics[@]} ; do
    echo "Producing input for systematic: "$systematic
    ./install/bin/Histo -d stacked -g ttbb -c combined -s ${systematic} 
done

echo "Production of all inputs finished"
echo
