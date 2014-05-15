#!/bin/sh

source $(dirname `readlink -f $0`)/parallelTools.sh

# in the excludeList, put all distributions that you dont want. Separate them with a |, i.e. HypLLBarDPhi|HypNeutrinopT
# excludeList='HypNeutrinopT|bcp|_step|akr|bkr|SF|#'
# plotList=`awk '{print $1}' HistoList_control | grep Hyp| grep -Erv $excludeList`

echo
echo "******************************************************************"
echo "Creating control plots for signal variations ONLY"
echo "If you want to create control plots for all/any other variation please run"
echo "install/bin/Histo -t cp -s <systematic>"
echo "******************************************************************"
echo
echo


excludeList="^$|step1|step2|step3|step4|akr|bkr|SF|#|bcp"
plotList=`awk '{print $1}' HistoList_control | grep -Ev $excludeList`


echo "Please press any key to start unfolding the following distributions in parallel or press Ctrl-C to cancel:"
echo "$plotList" | perl -l40 -pe ''
read -n 1 -s
echo ""

########################
## draw all systematic variations for control plots
## no unceratinty band will be plotted yet
########################
for i in $plotList; do
    $HISTO -t cp -s MATCH_UP -s MATCH_DOWN -s SCALE_UP -s SCALE_DOWN -s MCATNLO -s POWHEG -p +$i &
    $HISTO -t cp -s JES_UP -s JES_DOWN -s JER_UP -s JER_DOWN -p +$i &
    $HISTO -t cp -s PU_UP -s PU_DOWN -s LEPT_UP -s LEPT_DOWN -s TRIG_UP -s TRIG_DOWN -p +$i &
    $HISTO -t cp -s BTAG_UP -s BTAG_DOWN -s BTAG_PT_UP -s BTAG_PT_DOWN -s BTAG_ETA_UP -s BTAG_ETA_DOWN -p +$i &
    $HISTO -t cp -s BTAG_LJET_UP -s BTAG_LJET_DOWN -s BTAG_LJET_PT_UP -s BTAG_LJET_PT_DOWN -s BTAG_LJET_ETA_UP -s BTAG_LJET_ETA_DOWN -p +$i &
    w
done

########################
## draw Nominal control plot including uncertainty band
##   plese notice the ' -b'
########################

echo "----------------------------------------------------------------"
echo "Now sumbitting jobs to draw the control plot with the error band"
for i in $plotList; do 
    $HISTO -t cp -s Nominal -p +$i  -b&
    w
done
