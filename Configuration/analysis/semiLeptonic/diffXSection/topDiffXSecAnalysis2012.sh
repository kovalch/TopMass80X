#!/bin/sh

## this shell script runs all macro parts of the 2012 differential cross section analysis 

## HOW TO USE THIS SHELL SCRIPT
## note: has to be called from .../TopAnalysis/Configuration/analysis/semiLeptonic/diffXSection
## a) set up the folder structure:
## mkdir -p diffXSecFromSignal/plots/muon/2012/monitoring/withoutRatioPlots
## mkdir -p diffXSecFromSignal/plots/muon/2012/partonLevel
## mkdir -p diffXSecFromSignal/plots/muon/2012/recoYield
## mkdir -p diffXSecFromSignal/plots/muon/2012/uncertainties
## mkdir -p diffXSecFromSignal/plots/muon/2012/uncertaintyDistributions
## mkdir -p diffXSecFromSignal/plots/muon/2012/uncertaintyDistributionsOverview
## mkdir -p diffXSecFromSignal/plots/muon/2012/xSec
## mkdir -p diffXSecFromSignal/plots/muon/2012/binning
## mkdir -p diffXSecFromSignal/plots/muon/2012/effANdacc
## mkdir -p diffXSecFromSignal/plots/muon/2012/genRecoCorrPlots
## mkdir -p diffXSecFromSignal/plots/muon/2012/kinFitPerformance
## mkdir -p diffXSecFromSignal/plots/muon/2012/shapeReweighting
## mkdir -p diffXSecFromSignal/plots/muon/2012/unfolding
## mkdir -p diffXSecFromSignal/plots/muon/2012/closureTest
## mkdir -p diffXSecFromSignal/plots/muon/2012/topPtTest
## mkdir -p diffXSecFromSignal/plots/electron/2012/monitoring/withoutRatioPlots
## mkdir -p diffXSecFromSignal/plots/electron/2012/partonLevel
## mkdir -p diffXSecFromSignal/plots/electron/2012/recoYield
## mkdir -p diffXSecFromSignal/plots/electron/2012/uncertainties
## mkdir -p diffXSecFromSignal/plots/electron/2012/uncertaintyDistributions
## mkdir -p diffXSecFromSignal/plots/electron/2012/uncertaintyDistributionsOverview
## mkdir -p diffXSecFromSignal/plots/electron/2012/xSec
## mkdir -p diffXSecFromSignal/plots/electron/2012/binning
## mkdir -p diffXSecFromSignal/plots/electron/2012/effAndAcc
## mkdir -p diffXSecFromSignal/plots/electron/2012/genRecoCorrPlots
## mkdir -p diffXSecFromSignal/plots/electron/2012/kinFitPerformance
## mkdir -p diffXSecFromSignal/plots/electron/2012/shapeReweighting 
## mkdir -p diffXSecFromSignal/plots/electron/2012/unfolding
## mkdir -p diffXSecFromSignal/plots/electron/2012/closureTest
## mkdir -p diffXSecFromSignal/plots/electron/2012/topPtTest
## mkdir -p diffXSecFromSignal/plots/combined/2012/xSec
## mkdir -p diffXSecFromSignal/plots/combined/2012/binning
## mkdir -p diffXSecFromSignal/plots/combined/2012/partonLevel
## mkdir -p diffXSecFromSignal/plots/combined/2012/recoYield
## mkdir -p diffXSecFromSignal/plots/combined/2012/binning
## mkdir -p diffXSecFromSignal/plots/combined/2012/effAndAcc
## mkdir -p diffXSecFromSignal/plots/combined/2012/genRecoCorrPlots
## mkdir -p diffXSecFromSignal/plots/combined/2012/kinFitPerformance
## mkdir -p diffXSecFromSignal/plots/combined/2012/uncertainties
## mkdir -p diffXSecFromSignal/plots/combined/2012/uncertaintyDistributions
## mkdir -p diffXSecFromSignal/plots/combined/2012/uncertaintyDistributionsOverview
## mkdir -p diffXSecFromSignal/plots/combined/2012/monitoring/withoutRatioPlots
## mkdir -p diffXSecFromSignal/plots/combined/2012/xSecRatio
## mkdir -p diffXSecFromSignal/plots/combined/2012/effAndAcc/withoutRatioPlots
## mkdir -p diffXSecFromSignal/plots/combined/2012/ttgencomparison
## mkdir -p diffXSecFromSignal/plots/combined/2012/comparisonATLAS
## mkdir -p diffXSecFromSignal/plots/combined/2012/unfolding
## mkdir -p diffXSecFromSignal/plots/combined/2012/closureTest
## mkdir -p diffXSecFromSignal/plots/combined/2012/massConstraintTest
## mkdir -p diffXSecFromSignal/plots/combined/2012/regularizationTest
## mkdir -p diffXSecFromSignal/plots/combined/2012/topPtTest
## mkdir -p diffXSecFromSignal/plots/combined/2012/cov
## mkdir -p diffXSecFromSignal/plots/combined/2012/rhos

## b) root files needed for the Analysis are loaded automatically from /afs/naf.desy.de/group/cms/scratch/tophh/
## c) if not yet done, combine the MC samples for the single channels (like QCD, single top, Diboson) using combineMCsamples.C
## d) when using the shell script for the very first time, do "chmod a+x topDiffXSecAnalysis2012.sh
## find final plots in ./diffXSecFromSignal/plots/ after running the analysis via ./topDiffXSecAnalysis2012.sh
## e1) Specify the parameters under 'configure settings' in this shell script
## e2) The decay channel, PS and level are specified via 3 external arguments, e.g.
##    ./topDiffXSecAnalysis2012.sh muon (or electron or combined) extrapolate (or visible) parton (or hadron)
##    if 0 arguments are given -> default values as specified below
##    if not 0 or 3 arguments -> Abort!!!
## f) Further plots can be produced using the following macros:
##    ATLASCMSCOMPARISON.C -> 7TeV comparison plots CMS&ATLAS
##    btagEfficiencyComparison.C -> some b-tag related simulation comparison plots
##    getMCcorrectionfactors.C -> produce data/MadGraph scale factors from CMS data results
##    analyzeFileComparison.C -> compare mass bias from choice of top mass in reconstruction
##    applyBCCs.C -> produce BCC values and store them in a rootfile to be used in combineTopDiffXSecUncertainties
##    optimizeProbCut.C -> calculate optimal value for probability selection cut and produce some control ditsributions and numbers
##    createPseudoData.C -> create pseudo data based on the MC simluation (also mixing in e.g. Z')
##    addDistributions.C/addHistograms.C -> add additional control plots to the analyses root files
##    treeComparison.C-> plot topPt in bins of e.g. chi2 and nPV
##    ATLASCompTreeSGsamples.C-> compare various official ATLAS and CMS samples
##    analyzeGenComparison7TeV.C -> like analyzeGenComparison.C but for 7TeV samples/analysis
##    createNNLOplot.C -> create ahrens theory curves
##    analyzeRegularizationTest.C-> test different choices of regularization methods and parameters

########################
## configure settings ##
########################

## default values

# lepton flavour in semi leptonic decay
# choose \"muon\" or \"electron\" or \"combined\"
decayChannel=\"combined\" 

## extrapolate xSec to full PS?
## extrapolate = true / false (default: true)
extrapolate=true

## use hadron level PS instead of parton level PS?
## hadron = true / false (default: false)
hadron=false

## Unfolding closure test -> use pseudo data
## closureTestSpecifier = \"NoDistort\", \"topPtUp\", \"topPtDown\", \"ttbarMassUp\", \"ttbarMassDown\", \"data\" or \"1000\"
## for corresponding ttbarReweighting/Zprime pseudo data  (default: \"\")
closureTestSpecifier=\"\"

## run combination based on event yield combination instead of 
## combinedEventYields = true / false (default: false)
## attention: this affects only bothDecayChannelsCombination.C 
##            and is automatically adjusted with input string "combined2"
combinedEventYields=false

## use prob sel analysis
#addSel=\"\"
#chi2Max=99999
#probMin=-1
addSel=\"ProbSel\"
chi2Max=7.824 # corresponds to probMin
probMin=0.02  # corresponds to chi2Max
ndof=2 #corresponds to prob cut on second KinFit with equal mass constraint

## take arguments
clear
echo "-------------------------------------------------------------------------------------"
if [ $# -eq 0 ]; then
   echo "The default settings for decay channel, phase space extrapolation and level are used:"
elif [ $# -eq 3 -o $# -eq 4 ]; then
   ## set first argument
   if   [ $1 == "muon"      ];   then decayChannel=\"muon\"
   elif [ $1 == "electron"  ];   then decayChannel=\"electron\"
   elif [ $1 == "combined"  ];   then decayChannel=\"combined\"
   elif [ $1 == "combined2" ];   then 
       decayChannel=\"combined\"
       combinedEventYields=true
   else                          echo "1st argument ( $1 ) is not valid! Choose \"muon\", \"electron\" or \"combined\" as decay channel! Abort!"; exit
   fi
   ## set second argument
   if   [ $2 == "extrapolate" ]; then extrapolate=true
   elif [ $2 == "visible"     ]; then extrapolate=false
   else                          echo "2nd argument ( $2 ) is not valid! Choose \"extrapolate\" or \"visible\" as phase space! Abort!"; exit
   fi
   ## set third argument
   if   [ $3 == "parton" ];      then hadron=false
   elif [ $3 == "hadron" ];      then hadron=true
   else                          echo "3rd argument ( $3 ) is not valid! Choose \"parton\" or \"hadron\" as level! Abort!"; exit
   fi
## in case of fourth argument:
   if [ $# -eq 4 ]; then
       closureTestSpecifier=\"$4\"
   fi
   echo "Decay channel, phase space extrapolation and level were specified by user:"
else
   echo "Wrong number of arguments! Choose 0 arguments for default values and 3 arguments to specify" 
   echo "decay channel, phase space extrapolation and level. Abort!"
   exit
fi

echo "-------------------------------------------------------------------------------------"
echo "Decay channel:                              $decayChannel    "
echo "extrapolate:                                $extrapolate     "
echo "hadron:                                     $hadron          "
echo "closureTestSpecifier                        $closureTestSpecifier          "
echo "-------------------------------------------------------------------------------------"
echo

## folder on /afs/naf.desy.de/group/cms/scratch/tophh where MC and data files are stored
## inputFolderName=\"RecentAnalysisRun8TeV\" (default)
#inputFolderName=\"RecentAnalysisRun8TeV\"
inputFolderName=\"RecentAnalysisRun8TeV_doubleKinFit\"

## Dataset and luminosity [/pb]
## has to fit to current dataset

mudataSample=\"/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/muonDiffXSecData2012ABCDAll.root\"
eldataSample=\"/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/elecDiffXSecData2012ABCDAll.root\"

# closure tests
# [ $closureTestSpecifier != \"\" ]; then
#  if [ $closureTestSpecifierS==\"NoDistort\"]; then
#	mudataSample=\"/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/pseudodata/muonPseudoData19712pb8TeV.root\"
#	eldataSample=\"/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/pseudodata/electronPseudoData19712pb8TeV.root\"
#  elif [ $closureTestSpecifierS==\"topPtUp\"]; then
#	mudataSample=\"/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/pseudodata/muonPseudoData19712pbReweightedtopPtUp8TeV.root\"
#	eldataSample=\"/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/pseudodata/electronPseudoData19712pbReweightedtopPtUp8TeV.root\"
#  elif [ $closureTestSpecifierS==\"topPtDown\"]; then
#	mudataSample=\"/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/pseudodata/muonPseudoData19712pbReweightedtopPtDown8TeV.root\"
#	eldataSample=\"/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/pseudodata/electronPseudoData19712pbReweightedtopPtDown8TeV.root\"
#  elif [ $closureTestSpecifierS==\"ttbarMassUp\"]; then
#	mudataSample=\"/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/pseudodata/muonPseudoData19712pbReweightedttbarMassUp8TeV.root\"
#	eldataSample=\"/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/pseudodata/electronPseudoData19712pbReweightedttbarMassUp8TeV.root\"
#  elif [ $closureTestSpecifierS==\"ttbarMassDown\"]; then
#	mudataSample=\"/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/pseudodata/muonPseudoData19712pbReweightedttbarMassDown8TeV.root\"
#	eldataSample=\"/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/pseudodata/electronPseudoData19712pbReweightedttbarMassDown8TeV.root\"
#  elif [ $closureTestSpecifierS==\"data\"]; then
#	mudataSample=\"/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/pseudodata/muonPseudoData19712pbReweighteddata8TeV.root\"
#	eldataSample=\"/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/pseudodata/electronPseudoData19712pbReweighteddata8TeV.root\"
#  elif [ $closureTestSpecifierS==\"1000\"]; then
#	mudataSample=\"/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/pseudodata/muonPseudoData19712pbandM1000W100Zprime8TeV.root\"
#	eldataSample=\"/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/pseudodata/electronPseudoData19712pbandM1000W100Zprime8TeV.root\"
#  fi
#
# -> directly specified within basicFunctions.h

if [ $decayChannel == \"electron\" ]; then
    dataLuminosity=19712
    dataSample=$eldataSample
else
    if [ $decayChannel == \"muon\" ]; then
	dataLuminosity=19712
	dataSample=$mudataSample
    else
	dataLuminosity=19712 # mean value
	dataSample=$eldataSample\":\"$mudataSample
    fi
fi

## Data label, required for filename
## dataLabel=2012 (default)
dataLabel=2012

## save all plots? (.eps and .root)
save=true

## detail level of output 
## 0: no output (default)
## 1: std output
## 2: output for debugging
verbose=0

## Re-create monitoring plots
## redoControlPlots = true / false (default: true)
redoControlPlots=true

## ttbar generator comparison plots
## genComparison = true / false (default: true)
genComparison=true

## Re-do study for optimal value of the probability cut
## redoProb = true / false (default: true)
redoProb=true

## Re-create systematically varied results
## redoSystematics = true / false (default: true)
redoSystematics=true

## Re-create systematically PDF shifted signal files 
## redoPDFReweighting = true / false (default: true)
redoPDFReweighting=true

## Produce final xSec plots, ratios, uncertainties and covariance matrices 
## produceResults = true / false (default: true)
produceResults=true

## Produce covariance matrices for final result
## redoCov = true / false (default: true)
redoCov=true

## scan different regularization parameter
## regTest = true / false (default: true)
regTest=true

## scan different mass constraints
## constraintTest = true / false (default: true)
constraintTest=true

## check topPt in different slices
## topPtTests = true / false (default: true)
topPtTests=true

## Make pt plots logarithmic
## makeLogPlots = true / false (default: false)
makeLogPlots=false

## last systematic to proceed (0: only std analysis without variation)
## has to be consistent with the enumerator "systematicVariation" in "basicFunctions.h"
## maxSys>0 needs a lot of time (must be <= 59 (default), see list of systematics below)
maxSys=59

## Include cross-check variables to get additional differential cross-sections for
## a) pT(top) and pT(antitop)
## b) y(top)  and y(antitop)
## c) eta(l+) and eta(l-)
## Attention: The exectution mainly of analyzeHypothesisKinFit.C lasts longer if this parameter is set to true
##
## inclCCVars = true / false (default: false)
inclCCVars=false

## disable waiting time to read output
## fast = true / false (default: true)
fast=true

## delete all (old) existing .eps, .png and .pdf plots?
## clean = true / false (default: false)
clean=false

## use SVD unfolding?
## SVD = true / false (default: true)
SVD=true

## redetermine regularisation parameter tau?
## redetTau = true / false (default: false)
redetTau=false

## Re-create purity/stability/resolution plots
## redoPurStab = true / false (default: true)
redoPurStab=true

## Use bin-centre corrections (BCC)
## useBCC = true / false (default: false)
useBCC=false

#### =====================
####  Prepare running
#### =====================

BoolArray=( true false )

## start the timer to stop the time needed for running the whole analysis
START=$(date +%s)

# get proper name of rootfile
PS=""
LV="Parton"
if [ $extrapolate == false ] 
    then 
    PS="PhaseSpace"
    echo $PS
fi
if [ $extrapolate == false -a $hadron == true ] 
    then
    LV="Hadron"
fi

# switches for unfolding closure test
if [ $closureTestSpecifier != \"\" ] 
    then
    echo
    echo "CLOSURE TEST FOR UNFOLDING!"
    echo "closure test type: " $closureTestSpecifier
    redoPDFReweighting=false
    redoProb=false
    redoControlPlots=false
    redoSystematics=false
    maxSys=0
    redoPurStab=false
    redetTau=false
    genComparison=false
    useBCC=false
    clean=false
    redoCov=false
    topPtTests=false
    constraintTest=false
fi

muonFile=./diffXSecTopSemiMu$dataLabel$LV$PS.root
elecFile=./diffXSecTopSemiElec$dataLabel$LV$PS.root
combFile=./diffXSecTopSemiLep$dataLabel$LV$PS.root

## print out configuration
if [ $decayChannel == \"combined\" ]
    then
    echo
    echo "combining the electron and muon channel"
    if [ -f $muonFile -a -f $elecFile -a $1 != "combined2" ]; then
	echo
	echo "Doing the full differential top xSec analysis. "
	echo 
	echo "Data Label:                                 $dataLabel       "
	echo "Used data:                                  $dataSample      "
	echo "Decay channel:                              $decayChannel    "
	echo "Luminosity:                                 $dataLuminosity  " 
	echo "Re-do control plots:                        $redoControlPlots" 
	echo "Re-do systematic uncertainties:             $redoSystematics "
	echo "Number of considered systematics:           $maxSys          "
	echo "additional selection:                       $addSel (chi2<$chi2Max)"
	echo "Save plots:                                 $save            " 
	echo
    elif [  $1 == "combined2" ]; then
	echo
	echo "Doing the full differential top xSec analysis. "
	echo "COMBINING THE CHANNELS AT EVENT YIELD LEVEL"
	echo 
	echo "Data Label:                                 $dataLabel       "
	echo "Used data:                                  $dataSample      "
	echo "Decay channel:                              $decayChannel    "
	echo "Luminosity:                                 $dataLuminosity  "
	echo "Re-do control plots:                        $redoControlPlots" 
	echo "Re-do probability study:                    $redoProb        "
	echo "Re-do systematic uncertainties:             $redoSystematics "
	echo "Re-do covariance matrix::                   $redoCov         "
	echo "Perform regularization parameter scan test: $regTest         "
        echo "Perform top mass constraint bias test:      $constraintTest  "
        echo "Perform generator comparison test:          $genComparison   "
        echo "Perform topPt in different slices test:     $topPtTests      "
	echo "Produce final results:                      $produceResults  "
	echo "Number of considered systematics:           $maxSys          "
	echo "additional selection:                       $addSel (chi2<$chi2Max)"
	echo "Save plots:                                 $save            "
	echo
    else
	echo
	echo "NOTE: The combination requires two files"
	echo "a) "$muonFile
        echo "b) "$elecFile
	echo
	echo "Please get them by running the e/mu channel first"
	echo
	exit
    fi
fi

if [ $fast == false ]
    then
    sleep 5
fi

#### ===================================
####  Delete existing root file/ plots
#### ===================================
echo
echo "part A: Delete existing files and plots (if applicable)"
if [ $fast == false ]; then
    sleep 3
fi

if [ $clean == true ]; then

    if [ $redoSystematics == false ]; then
	echo
	echo "Flag 'redoSystematics' is set to $redoSystematics "
	echo "Flag 'clean' set to 'false' to avoid deleting files/plots which are not recreated "
	echo
	clean=false
    else
	
        ## ============================
        ##  Delete existing root files
        ## ============================
   
	echo
	echo "Part A1: delete existing root file"
	if [ $decayChannel == \"electron\" ]; then
	    echo $elecFile
	    rm $elecFile
	else
	    if [ $decayChannel == \"muon\" ]; then
		echo $muonFile
		rm $muonFile
	    else
		if [ $decayChannel == \"combined\" ]; then
		    echo $muonFile
		    rm $combFile		    
		fi
	    fi
	fi 
	
        ## ============================    
        ##  Delete existing plots
	## ============================
    
	echo
	echo "Part A2: delete existing plots within diffXSecFromSignal/plots/$decayChannel/$dataLabel/*/*.*"
	if [ $decayChannel == \"electron\" ]; then
	    rm ./diffXSecFromSignal/plots/electron/$dataLabel/*/*.*
	else
	    if [ $decayChannel == \"muon\" ]; then
		rm ./diffXSecFromSignal/plots/muon/$dataLabel/*/*.*
	    else
		if [ $decayChannel == \"combined\" ]; then
		    rm ./diffXSecFromSignal/plots/combined/$dataLabel/*/*.*
		fi
	    fi
	fi
    fi
fi

#### =====================
####  Run cut monitoring
#### =====================

BEFOREB=$(date +%s)
echo
echo "Part B: process cut monitoring macros"
if [ $fast == false ]
    then
    sleep 3
fi

## prob cut
if [ $redoProb = true -a $redoControlPlots = true ]; then
    root -l -q -b './optimizeProbCut.C++g('$probMin', '$ndof', '$verbose', '$decayChannel', '$inputFolderName', '$save')'
fi

if [ $redoControlPlots == true ]; then
    
    ## Compile library
    
    if [ -f commandsMonPrepare.cint ]; then    
	rm commandsMonPrepare.cint
	rm analyzeTopDiffXSecMonitoring_C.so
	rm analyzeTopDiffXSecMonitoring_C.d
    fi
    
    cat >> commandsMonPrepare.cint << EOF
.L analyzeTopDiffXSecMonitoring.C++
EOF
    
    root -l -b < commandsMonPrepare.cint
    
    ## Execute macro

    if [ $redoControlPlots == true ]; then
	
	for label in "${BoolArray[@]}"; do
	    
        ## label = 1: Control plots with ratio plots
        ## label = 0: Control plots without ratio plots
	    
	    if [ -f commandsMonRun.cint ]; then    
		rm commandsMonRun.cint       
	    fi    
	    
	    cat >> commandsMonRun.cint << EOF
.L analyzeTopDiffXSecMonitoring_C.so
analyzeTopDiffXSecMonitoring($dataLuminosity, $save, $verbose, $inputFolderName, $dataSample, $decayChannel, $label, $extrapolate, $hadron, $addSel) 
EOF
	    echo ""
	    echo " Processing .... analyzeTopDiffXSecMonitoring.C++($dataLuminosity, $save, $verbose, $inputFolderName, $dataSample, $decayChannel, true, $extrapolate, $hadron)"
	    root -l -b < commandsMonRun.cint
	done
	echo ""
	echo " Processing .... analyzeGenComparison.C++($save', '$verbose, true/false)"
        ## gen comparison
	if [ $genComparison == true ]; then  
	    root -l -q -b './analyzeGenComparison.C++('$save', '$verbose', true )'
	    root -l -q -b './analyzeGenComparison.C++('$save', '$verbose', false)'	
	fi
    fi
fi

### ===================================
###  Run migration macro for binning 
### ===================================

BEFOREC=$(date +%s)
echo
echo "Part C: process migration macro to validate binning"
if [ $fast == false ]
    then
    sleep 3
fi

if [ $redoPurStab == true ]
    then
    # Array of differential variables
    listVar_=( \"topPt\" \"topY\" \"ttbarPt\" \"ttbarY\" \"ttbarMass\" \"topPtTtbarSys\" \"ttbarDelPhi\" \"ttbarPhiStar\" )
    plotAcceptance=true
    if [ $hadron == true ]; 
	then 
	listVar_=( \"lepPt\" \"lepEta\" \"bqPt\" \"bqEta\" \"bbbarPt\" \"bbbarMass\" \"lbMass\" \"rhos\" \"Njets\" )
	Plotacceptance=false
    fi
    	
    echo "purity and stability will be calculated for the following variables: "
    echo
    echo "${listVar_[@]}"
  
    # loop over all systematic variations
    for (( iVar=0; iVar<${#listVar_[@]}; iVar++ )); do
	root -l -q -b './purityStabilityEfficiency.C++('${listVar_[$iVar]}','$save', '$decayChannel', '$inputFolderName', '$plotAcceptance', true, false, '$chi2Max', 1, '$hadron')'
    done
fi



#### ============================
####  Prepare PDF uncertainties 
#### ============================
BEFORED=$(date +%s)
echo
echo "Part D: Prepare files for pdf uncertainties"

if [ $decayChannel != \"combined\" -a $redoSystematics == true -a $redoPDFReweighting == true ]; then
    echo
    root -l -q -b './analyzeTopDiffXSecMCdependency.C++('$dataLuminosity','$decayChannel', '$save', '$verbose', '$inputFolderName', '$dataSample', 'true', '$inclCCVars')' 
elif [ $1 == "combined2" -a $redoSystematics == true -a $redoPDFReweighting == true ]; then
    root -l -q -b './analyzeTopDiffXSecMCdependency.C++('$dataLuminosity', '\"muon\"',     '$save', '$verbose', '$inputFolderName', '$mudataSample', 'true', '$inclCCVars')' 
    root -l -q -b './analyzeTopDiffXSecMCdependency.C++('$dataLuminosity', '\"electron\"', '$save', '$verbose', '$inputFolderName', '$eldataSample', 'true', '$inclCCVars')' 
else
    echo "Done for 2012 analysis (in e/mu channel separate or when combining event yields) and if systematics and PDF files are requested to be re-done (redoSystematics set to $redoSystematics, redoPDFReweighting set to $redoPDFReweighting)."
fi

#### ==========================================
####  Run efficiency & cross section macro 
####  ==========================================
BEFOREE=$(date +%s)
echo
echo "part E1: process cross section calculation macro for all systematics"
echo "INFO: missing files must not be problematic"
echo "      either all WZ, WW and ZZ or the combined VV sample are necessary"
echo "      same is true for the single top samples (s, t, tW)"
echo
if [ $fast == false ]
    then
    sleep 5
fi

## print key for systematic variations
## has to be consistent with the enumerator "systematicVariation" in "basicFunctions.h"

echo

echo "  0: sysNo                                                   "                                                    
echo "  1: sysLumiUp                   2: sysLumiDown              "  
echo "  3: sysPUUp                     4: sysPUDown                "  
echo "  5: sysJESUp                    6: sysJESDown               "  
echo "  7: sysJERUp                    8: sysJERDown               "  
echo "  9: sysLepEffSFNormUp          10: sysLepEffSFNormDown      "  
echo " 11: sysLepEffSFShapeUpEta      12: sysLepEffSFShapeDownEta  "  
echo " 13: sysLepEffSFShapeUpPt       14: sysLepEffSFShapeDownPt   "  
echo " 15: sysBtagSFUp                16: sysBtagSFDown            "  
echo " 17: sysBtagSFShapeUpPt65       18: sysBtagSFShapeDownPt65   "  
echo " 19: sysBtagSFShapeUpEta0p7     20: sysBtagSFShapeDownEta0p7 "  
echo " 21: sysMisTagSFUp              22: sysMisTagSFDown          "  
echo " 23: sysTopScaleUp              24: sysTopScaleDown          "  
echo " 25: sysVBosonScaleUp           26: sysVBosonScaleDown       "  
echo " 27: sysSingleTopScaleUp        28: sysSingleTopScaleDown    "  
echo " 29: sysTopMatchUp              30: sysTopMatchDown          "  
echo " 31: sysVBosonMatchUp           32: sysVBosonMatchDown  	   "
echo " 33: sysTopMassUp               34: sysTopMassDown  	   "
echo " 35: sysTopMassUp2              36: sysTopMassDown2	   "
echo " 37: sysTopMassUp3              38: sysTopMassDown3	   "
echo " 39: sysTopMassUp4              40: sysTopMassDown4	   "
echo " 41: sysQCDUp                   42: sysQCDDown               "  
echo " 43: sysSTopUp                  44: sysSTopDown              "  
echo " 45: sysDiBosUp                 46: sysDiBosDown 		   "
echo " 47: sysVjetsUp                 48: sysVjetsDown		   "
echo " 49: sysBRUp                    50: sysBRDown                "
echo " 51: sysPDFUp                   52: sysPDFDown               "  
echo " 53: sysUnf                     54: sysMad                   "
echo " 55: sysHadUp                   56: sysHadDown               "
echo " 57: sysGenMCatNLO              58: sysGenPowheg  	   "
echo " 59: sysGenPowhegHerwig         60: ENDOFSYSENUM             "

echo

if [ $fast == false ]; then
    sleep 5
fi

#### ============================
####  Systematic Uncertainties
#### ============================

BEFORESYS=$(date +%s)

if [ $decayChannel != \"combined\" -o $1 == "combined2" ]; then

    ## ====================================================================================
    ##  Compile library, required only once before processing systematic uncertainties
    ## ====================================================================================
    
    if [ -f commandsSysPrepare.cint ]; then    
	rm commandsSysPrepare.cint
	rm analyzeHypothesisKinFit_C.so
	rm analyzeHypothesisKinFit_C.d
    fi

    cat >> commandsSysPrepare.cint << EOF
.L analyzeHypothesisKinFit.C++g
EOF

    root -l -b < commandsSysPrepare.cint

    ## ================================================================================================================
    ##  Processing reference data (noSys), required always thus excluded from looping over systematic uncertainties 
    ## ================================================================================================================

    if [ -f commandsNoSysRun.cint ]; then    
	rm commandsNoSysRun.cint
    fi
    
    cat >> commandsNoSysRun.cint << EOF
.L analyzeHypothesisKinFit_C.so
analyzeHypothesisKinFit($dataLuminosity, $save, 0, $verbose, $inputFolderName, $dataSample, $decayChannel, $SVD, $extrapolate, $hadron, $inclCCVars, $redetTau, $closureTestSpecifier, $addSel)
EOF

    echo ""
    echo " Processing .... analyzeHypothesisKinFit($dataLuminosity, $save, 0, $verbose, $inputFolderName, $dataSample, $decayChannel, $SVD, $extrapolate, $hadron ,$inclCCVars, $redetTau, $closureTestSpecifier, $addSel)"
    root -l -b < commandsNoSysRun.cint


    ## ==========================================
    ##  Processing systematic uncertainties 
    ## ==========================================

    if [ $redoSystematics == true ]; then
    
        ## loop all systematic variations
	
	for (( systematicVariation = 1; systematicVariation <= $maxSys;  systematicVariation++ )); do

	    ## run macro for 2012 analysis
	    
	    if [ -f commandsSysRun.cint ]; then    
		rm commandsSysRun.cint
	    fi
	    
	    cat >> commandsSysRun.cint << EOF
.L analyzeHypothesisKinFit_C.so
analyzeHypothesisKinFit($dataLuminosity, $save, $systematicVariation, $verbose, $inputFolderName, $dataSample, $decayChannel, $SVD, $extrapolate, $hadron, $inclCCVars, $redetTau, $closureTestSpecifier, $addSel)
EOF
	    echo ""
	    echo " Processing .... analyzeHypothesisKinFit($dataLuminosity, $save, $systematicVariation, $verbose, $inputFolderName, $dataSample, $decayChannel, $SVD, $extrapolate, $hadron, $inclCCVars, $redetTau, $closureTestSpecifier, $addSel)"
	    root -l -b < commandsSysRun.cint
	done
    fi
else
    echo "will be ignored, only done for decayChannel=muon/electron"
fi

#### ===================================
####  Combine electron and muon channel 
#### ===================================
echo
echo "Part E2: Combine electron and muon channel"
if [ $fast == false ]; then
    sleep 2
fi

if [ $decayChannel == \"combined\" -a $produceResults == true ]; then
    
    ## A final uncertainties and result plots
    echo "Cross sections for all systematic variations and combined decay channels"
    echo
    ## delete old files    
    if [ -f commandsCombineChannelsPrepare.cint ]; then    
	rm commandsCombineChannelsPrepare.cint
	rm bothDecayChannelsCombination_C.so
	rm bothDecayChannelsCombination_C.d
    fi

    ## Compile library
    cat >> commandsCombineChannelsPrepare.cint << EOF
.L bothDecayChannelsCombination.C++g
EOF
    root -l -b < commandsCombineChannelsPrepare.cint
    
    ## Execute macros
    
    if [ -f commandsCombineChannelsRun.cint ]; then    
	rm commandsCombineChannelsRun.cint       
    fi  
    cat >> commandsCombineChannelsRun.cint << EOF
.L bothDecayChannelsCombination_C.so
bothDecayChannelsCombination($dataLuminosity, $save, $verbose, $inputFolderName, $makeLogPlots, $extrapolate, $hadron, $inclCCVars, $combinedEventYields, $closureTestSpecifier, $useBCC)
EOF

    echo ""
    echo " Processing .... bothDecayChannelsCombination($dataLuminosity, $save, $verbose, $inputFolderName, $makeLogPlots, $extrapolate, $hadron, $inclCCVars, $combinedEventYields, $closureTestSpecifier, $useBCC)"    
    root -l -b < commandsCombineChannelsRun.cint

    ## B covariance matrix
    echo
    echo "Covariance Matrix for systematic variations"
    echo
    if [ $redoCov == true ]; then    
        ## delete old files
	if [ -f commandsCovMatrixPrepare.cint ]; then    
	rm commandsCovMatrixPrepare.cint
	rm covarianceOfSystematicUnc_C.d
	rm covarianceOfSystematicUnc_C.so
	fi
        ## Compile library
	cat >> commandsCovMatrixPrepare.cint << EOF
.L covarianceOfSystematicUnc.C++g
EOF
	root -l -b < commandsCovMatrixPrepare.cint
	
        ## Execute macros    
	if [ -f commandsCovMatrixRun.cint ]; then    
	    rm commandsCovMatrixRun.cint       
	fi    
    
	cat >> commandsCovMatrixRun.cint << EOF
.L  covarianceOfSystematicUnc_C.so
covarianceOfSystematicUnc($save, $verbose+1, $decayChannel, $extrapolate, $hadron, $closureTestSpecifier)
EOF
	echo ""
	echo " Processing ....  covarianceOfSystematicUnc($save, $verbose, $decayChannel, $extrapolate, $hadron, $closureTestSpecifier)"
	root -l -b < commandsCovMatrixRun.cint
    else
	echo "... skipped"
    fi
else
    echo "will be ignored, only done for decayChannel=combined and if final results are produced (produceResults is set to $produceResults)"
fi
echo

AFTERSYS=$(date +%s)

#### ==========================================
####  Combine uncertainties for final xSecs 
#### ==========================================

BEFOREF=$(date +%s)

echo
echo "Part F: Calculate systematic errors and draw final cross section"
if [ $fast == false ]; then
    sleep 3
fi

if [ $produceResults == true -a $closureTestSpecifier == \"\" ]; then
    echo ""
    echo " Processing .... combineTopDiffXSecUncertainties($dataLuminosity, $save, $verbose, $decayChannel, $extrapolate, $hadron, $inclCCVars, $useBCC)"
    root -l -q -b './combineTopDiffXSecUncertainties.C++('$dataLuminosity', '$save', '$verbose', '$decayChannel', '$extrapolate', '$hadron', '$inclCCVars', '$closureTestSpecifier', '$useBCC')'
else
    echo "will be ignored, only done if final results are produced (produceResults is set to $produceResults)"
fi  

#### ==========================================
####  Create ratio plots for final xSecs 
#### ==========================================

if [ $decayChannel == \"combined\" -a $closureTestSpecifier == \"\" -a $produceResults == true ]; then
    echo ""
    echo " Processing .... createTheoryDataRatios($extrapolate, $hadron, $verbose)"
    root -l -q -b './createTheoryDataRatios.C++('$extrapolate', '$hadron', '$verbose')'
else
    echo "will be ignored, only done if final results are produced (produceResults is set to $produceResults)"
fi


#### ===================================================
####  Create latex code result tables for final xSecs 
#### ===================================================
if [ $decayChannel == \"combined\" -a $closureTestSpecifier == \"\" -a $produceResults == true ]; then
    echo ""
    echo " Processing .... makeResultTables($decayChannel, $extrapolate, $hadron, $inclCCVars, $useBCC)"
    root -l -q -b './makeResultTables.C++('$decayChannel', '$extrapolate', '$hadron', '$inclCCVars, $useBCC)'
else
    echo "will be ignored, only done if final results are produced (produceResults is set to $produceResults)"
fi


#### ===================================================
####  Perform regularization parameter scan
#### ===================================================
BEFOREG=$(date +%s)

echo
echo "Part G: Perform different studies and tests"

## collect tests to be performed 
listRegScan_=( \"mix\" )
if [ $closureTestSpecifier == \"\" ];
    then
    listRegScan_=( \"mix\" \"tau\" \"k\" )
fi
echo ""
echo " Processing .... analyzeRegularizationTest.C++(i, $dataLuminosity, $save, 0, $verbose, $inputFolderName, $dataSample, $decayChannel, $SVD, $extrapolate, $hadron, $redetTau, $closureTestSpecifier, $addSel)"
if [ $regTest == true -a $hadron == false -a $extrapolate == true ]; then 
    # loop over all tests
    for (( iTest=0; iTest<${#listRegScan_[@]}; iTest++ )); do
	echo "  - i=${listRegScan_[$iTest]}"
	root -l -b -q './analyzeRegularizationTest.C++g('${listRegScan_[$iTest]}','$dataLuminosity', '$save', 0, '$verbose', '$inputFolderName', '$dataSample', '$decayChannel', '$SVD', '$extrapolate', '$hadron', '$redetTau', '$closureTestSpecifier', '$addSel')'
    done
else
     echo "will be ignored, only done for extrapolate=true, hadron=false and regTest=true"
fi


#### ===================================================
####  Perform mass parameter bias test
#### ===================================================
echo ""
echo " Processing .... analyzeFileComparison($save, true/false , $verbose)"

if [ $constraintTest == true -a $closureTestSpecifier == \"\" ]; then
    root -l -b -q './analyzeFileComparison.C++g('$save', true , '$verbose')'
    root -l -b -q './analyzeFileComparison.C++g('$save', false, '$verbose')'
else
    echo "will be ignored (only done for real data and constraintTest=true)"
fi

#### ===================================================
####  Perform topPt in slices test
#### ===================================================
echo ""
echo " Processing .... treeComparison($dataLuminosity, $save, $verbose, $inputFolderName, $dataSample, $decayChannel, true, prob/PV)"

 listChecks_=( \"prob\" \"PV\" )
if [ $topPtTests == true -a $closureTestSpecifier == \"\" ]; then
    # loop over all tests
    for (( iCheck=0; iCheck<${#listChecks_[@]}; iCheck++ )); do
	echo "  - i=${listChecks_[$iCheck]}"
	root -l -b -q './treeComparison.C++g('$dataLuminosity', '$save', '$verbose', '$inputFolderName', '$dataSample', '$decayChannel', true, '${listChecks_[$iCheck]}')'
    done
else
    echo "will be ignored, only done for real data and topPtTests=true"
fi

#### ==========================================

echo ""
echo " All analysis steps finished!"
echo ""

#### ==========================================

#### ============================
####  After running the analysis
#### ============================
## stop the timer and echo time
END=$(date +%s)
TIME=$(( $END - $START ))
echo "time needed: $TIME seconds"
if [ $maxSys -ge 1 ]
    then
    SYS=$(( $AFTERSYS - $BEFORESYS ))
    echo "($SYS seconds due to systematic variations)"
fi
echo "part A: $(( $BEFOREB   - $START    )) seconds (clean up  )"
echo "part B: $(( $BEFOREC   - $BEFOREB  )) seconds (monitoring)"
echo "part C: $(( $BEFORED   - $BEFOREC  )) seconds (migration)"
echo "part D: $(( $BEFOREE   - $BEFORED  )) seconds (prepare PDF uncertainty run)"
echo "part F: $(( $BEFOREF   - $BEFOREE  )) seconds (xSecs- $maxSys systematic variations considered, cov matrix)"
echo "part G: $(( $BEFOREG   - $BEFOREF  )) seconds (errors and final xSec)"
echo "part H: $(( $END       - $BEFOREG  )) seconds (regularization parameter scan, mass constraint test"
