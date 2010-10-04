#!/bin/sh

## to be called from TopAnalysis/Configuration/analysis/semiLeptonic/diffXSection/

## use this file to do the whole analysis for:
## 1) l+jets cross section measurement [Njets]
## 2) inclusive top cross section measurement

## define processed data to be analyzed (output of data: ./analyzeMuonDiffXSec_cfg.py - MC: ./analyzeMuonDiffXSec_cfg.py)
dataSample=\"./diffXSecFromSignal/data/data0309/analyzeDiffXData_2900nb_residualJC.root\"
dataLuminosity=2880
echo
echo doing the full l+jets analysis
echo used data: $dataSample
echo luminosity: $dataLuminosity
echo
sleep 2
## delete old data plots
rm -r ./diffXSecFromSignal/plots/earlyData/*/*.*

## a) Cutflow and N-1 plots
## create plots for the analyzeMuonCuts.C Makro

echo processing analyzeMuonCuts.C
echo a\) create cutlow plots
sleep 2
date >> 'diffXSecFromSignal/plots/earlyData/cutMonitoring/cutflowNumbers.txt'
root -l -q -b './analyzeMuonCuts.C+('$dataLuminosity', true, '$dataSample', "cutflow")' >> 'diffXSecFromSignal/plots/earlyData/cutMonitoring/cutflowNumbers.txt'
echo b\) create N-1 plots
sleep 2
date >> 'diffXSecFromSignal/plots/earlyData/NminusOneDistributions/NminusOneNumbers.txt'
root -l -q -b './analyzeMuonCuts.C+('$dataLuminosity', true, '$dataSample', "NminusOne")' >> 'diffXSecFromSignal/plots/earlyData/NminusOneDistributions/NminusOneNumbers.txt'
echo analyzeMuonCuts.C ready

## b) do the bkg estimation and calculate MC effiencies 
##    -> all important numbers will be saved in crossSectionCalculation.txt
## example: makro.C+(luminosity, savePlots, writeInto:crossSectionCalculation.txt, pathOfDataFile.root)

rm ./crossSectionCalculation.txt
date >> './crossSectionCalculation.txt'
## b0) create pseudo data (needed in some makros)
echo create pseudo data for 50\/pb \(createPseudoData.C\)
sleep 2
root -l -q -b './createPseudoData.C+(50.)'
echo done
## b1) calculate R for charge asymmetry method
rm ./diffXSecFromSignal/plots/chargeAsymmetrie/Rcalculation.txt
echo processing chargeAsymmetrieCalculator.C
sleep 2
root -l -q -b './chargeAsymmetrieCalculator.C+(true, true)' >> './diffXSecFromSignal/plots/chargeAsymmetrie/Rcalculation.txt'
## b2) estimate N(QCD) via ABCD (MC - plots: 50 / pb)
rm ./diffXSecFromSignal/plots/ABCD/ABCDnumbers.txt
echo processing analyzeMuonDiffXABCDbkg.C \(plots done for 50\/pb\)
sleep 2
echo done
root -l -q -b './analyzeMuonDiffXABCDbkg.C+(50., true, true, '$dataSample')' >> './diffXSecFromSignal/plots/ABCD/ABCDnumbers.txt'
## b3) estimate MC based effiencies for l+jets correction (MC - plots: 5 / pb)
rm ./diffXSecFromSignal/plots/efficiency/efficiencyNumbers.txt
echo processing analyzeMuonDiffXEfficiency.C \(plots done for 5\/pb\)
sleep 2
root -l -q -b './analyzeMuonDiffXEfficiency.C+(5., true, true, false, "")' >> './diffXSecFromSignal/plots/efficiency/efficiencyNumbers.txt'
echo done
## b4) estimate N(W) via charge asymmetry method (MC - plots: 50 / pb)
rm ./diffXSecFromSignal/plots/chargeAsymmetrie/wjetsEstimationNumbers.txt
echo processing wjetsAsymmetrieEstimator.C \(plots done for 50\/pb\)
sleep 2
root -l -q -b './wjetsAsymmetrieEstimator.C+(50., true, true, '$dataSample')' >> './diffXSecFromSignal/plots/chargeAsymmetrie/wjetsEstimationNumbers.txt'
echo done

## c) systematic variations
## do the whole analysis applying different systematic variations
## and save all results in .txt file needed to calculate systematic errors
## example: analyzeMuonDiffXSec.C+(luminosity, savePlots, applyCorrections:.txtFile, pathOfDataFile.root, useMG=false, JES="", lumiShift=1.0, EffScaleFactor=1.0, QCDVariation=1.0, finalPlots=false)
echo processing systematic variations
echo 1 as reference without variation
sleep 2
root -l -q -b './analyzeMuonDiffXSec.C+('$dataLuminosity', true, true, '$dataSample', false, "", 1.0, 1.0, 1.0, 1.0, false)'
echo 2 using Madgraph ttbar absorped via eff
sleep 2
root -l -q -b './analyzeMuonDiffXEfficiency.C+(5., true, true, true, "")' >> './diffXSecFromSignal/plots/efficiency/efficiencyNumbers.txt'
root -l -q -b './analyzeMuonDiffXSec.C+('$dataLuminosity', true, true, '$dataSample', true, "", 1.0, 1.0, 1.0, 1.0, false)'
echo 3 JES shift up absorped via eff
sleep 2
root -l -q -b './analyzeMuonDiffXEfficiency.C+(5., true, true, false, "up")' >> './diffXSecFromSignal/plots/efficiency/efficiencyNumbers.txt'
root -l -q -b './analyzeMuonDiffXSec.C+('$dataLuminosity', true, true, '$dataSample', false, "up", 1.0, 1.0, 1.0, 1.0, false)'
echo 4 JES shift down absorped via eff
sleep 2
root -l -q -b './analyzeMuonDiffXEfficiency.C+(5., true, true, false, "down")' >> './diffXSecFromSignal/plots/efficiency/efficiencyNumbers.txt'
root -l -q -b './analyzeMuonDiffXSec.C+('$dataLuminosity', true, true, '$dataSample', false, "down", 1.0, 1.0, 1.0, 1.0, false)'
echo 5 Lumi shift up
sleep 2
root -l -q -b './analyzeMuonDiffXSec.C+('$dataLuminosity', true, true, '$dataSample', false, "", 1.11, 1.0, 1.0, 1.0, false)'
echo 6 Lumi shift down
sleep 2
root -l -q -b './analyzeMuonDiffXSec.C+('$dataLuminosity', true, true, '$dataSample', false, "", 0.89, 1.0, 1.0, 1.0, false)'
echo 7 Eff shift up
sleep 2
root -l -q -b './analyzeMuonDiffXSec.C+('$dataLuminosity', true, true, '$dataSample', false, "", 1.0, 1.045, 1.0, 1.0, false)'
echo 8 Eff shift down
sleep 2
root -l -q -b './analyzeMuonDiffXSec.C+('$dataLuminosity', true, true, '$dataSample', false, "", 1.0, 0.955, 1.0, 1.0, false)'
echo 9 QCD estimation up
sleep 2
root -l -q -b './analyzeMuonDiffXSec.C+('$dataLuminosity', true, true, '$dataSample', false, "", 1.0, 1.0, 2.0, 1.0, false)'
echo 10 QCD estimation down
sleep 2
root -l -q -b './analyzeMuonDiffXSec.C+('$dataLuminosity', true, true, '$dataSample', false, "", 1.0, 1.0, 0.0, 1.0, false)'
echo 11 Wjets estimation up
sleep 2
root -l -q -b './analyzeMuonDiffXSec.C+('$dataLuminosity', true, true, '$dataSample', false, "", 1.0, 1.0, 1.0, 1.3, false)'
echo 12 Wjets estimation down
sleep 2
root -l -q -b './analyzeMuonDiffXSec.C+('$dataLuminosity', true, true, '$dataSample', false, "", 1.0, 1.0, 1.0, 0.7, false)'

## d) do the final cross section calculation
## create final plots using the analyzeMuonDiffXSec.C Makro 
## calculate systematic errors and save all output within one single .txt file
## example: analyzeMuonDiffXSec.C+(luminosity, savePlots, applyCorrections:.txtFile, pathOfDataFile.root)
echo processing analyzeMuonDiffXSec.C
sleep 2
date >> 'diffXSecFromSignal/plots/earlyData/crossSection/diffXSecNumbers.txt'
root -l -q -b './analyzeMuonDiffXSec.C+('$dataLuminosity', true, true, '$dataSample', false, "", 1.0, 1.0, 1.0, 1.0, true, false)' | tee -a 'diffXSecFromSignal/plots/earlyData/crossSection/diffXSecNumbers.txt'
## with logarithmic plots
root -l -q -b './analyzeMuonDiffXSec.C+('$dataLuminosity', true, true, '$dataSample', false, "", 1.0, 1.0, 1.0, 1.0, true, true)' | tee -a 'diffXSecFromSignal/plots/earlyData/crossSection/diffXSecNumbers.txt'
echo

## e) monitor systematic shifts(JES, N(W), QCD(MC&ABCD) )
rm ./diffXSecFromSignal/plots/systematicVariations/*.*
echo monitor systematic shifts
sleep 2
## 1) in linear scale 
root -l -q -b './systematicUncertaintyScaling.C+('$dataLuminosity', true, '$dataSample', false)'
## 2) in logarithmic scale 
root -l -q -b './systematicUncertaintyScaling.C+('$dataLuminosity', true, '$dataSample', true)'

## f) copy plot to local folder on AFS
echo starting to copy plots...
## i)   analysis
scp -r ./diffXSecFromSignal/plots/earlyData/* mgoerner@uhh-cms014.desy.de:/afs/desy.de/user/m/mgoerner/private/note/pictures/
## ii)  ABCD results
scp -r ./diffXSecFromSignal/plots/ABCD mgoerner@uhh-cms014.desy.de:/afs/desy.de/user/m/mgoerner/private/note/pictures/
## iii) w+jets
scp -r ./diffXSecFromSignal/plots/chargeAsymmetrie mgoerner@uhh-cms014.desy.de:/afs/desy.de/user/m/mgoerner/private/note/pictures/
## iv)  efficiency
scp -r ./diffXSecFromSignal/plots/efficiency mgoerner@uhh-cms014.desy.de:/afs/desy.de/user/m/mgoerner/private/note/pictures/
## v)  efficiency
scp -r ./diffXSecFromSignal/plots/systematicVariations mgoerner@uhh-cms014.desy.de:/afs/desy.de/user/m/mgoerner/private/note/pictures/
echo ...ready