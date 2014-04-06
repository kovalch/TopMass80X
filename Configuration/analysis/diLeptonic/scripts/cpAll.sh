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

plotList="BjetEtaakr  BjetEtabkr   DIMFull_fullSel  DIMFull basic_jet_pt_step7 basic_lepton_pt_step7
          ElectronEta  ElectronEta_postMETcut  ElectronpT  ElectronpT_postMETcut
          HypBBBarMass  HypBJetEta  HypBjetMulti_noBTag  HypBjetMulti  HypBJetpT  HypjetMulti_diLep  HypjetMulti_noBTag  HypjetMulti  HypjetMultiXSec  HypLeptonBjetMass
          HypLeptonEta  HypLeptonpT  HypLLBarMass  HypLLBarpT  HypMet  HypTopMass  HypToppT  HypToppTTTRestFrame  HypTopRapidity  HypTTBarDeltaPhi  HypTTBarMass  HypTTBarpT  HypTTBarRapidity
          jetHT  jetpT  LeptonEtaakr  LeptonEtabkr  LeptonEta_diLep  LeptonEta  LeptonEta_postMETcut  LeptonpTakr  LeptonpTbkr  LeptonpT  LeptonpT_postMETcut
          METakr  METbkr  MET  MuonEta  MuonEta_postMETcut  MuonpT  MuonpT_postMETcut  step8  step9  triggersfeta  vertMulti_noPU  vertMulti"

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



########################
## draw Nominal control plot including uncertainty band
##   plese notice the ' -b'
########################
for i in $plotList; do 
    $HISTO -t cp -s Nominal -p +$i  -b&
    w
done
