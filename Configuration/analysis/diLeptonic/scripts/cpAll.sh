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
    $HISTO -t cp -s MATCH_UP -s MATCH_DOWN -s Nominal -s SCALE_UP -s SCALE_DOWN -s MCATNLO -s POWHEG -p +$i &
    w
done

## wait until all jobs are ready: wait for 5minutes
echo "Sleep 5mins until all control plots of the signal systematics are done"
sleep 5m

########################
## draw Nominal control plot including uncertainty band
##   plese notice the ' -b'
########################
for i in $plotList; do 
    $HISTO -t cp -s Nominal -p +$i  -b&
    w
done
