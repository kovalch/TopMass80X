#!/bin/bash

# This script is using the official Higgs analysis combine tool
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideHiggsAnalysisCombinedLimit

# USAGE: fitHfFrac scripts/fitHF.sh <folderName>

location=$PWD
cmssw_original=$CMSSW_BASE/src/
cmssw_combine="/data/group/top/HiggsCombineTool/CMSSW_6_1_1/src/"

if [ ! -d $cmssw_combine ]; then
    echo "### ERROR! No installed release with Higgs combine tool found at:"
    echo "###        $cmssw_combine"
fi

# Setting initial values for the input
workingFolder=$location"/HfFracScaleFactors"
if [ "$1" != "" ]; then
    workingFolder=$location"/$1"
fi
if [ ! -d $workingFolder  ]; then
    echo "### ERROR! Working directory doesn't exist: $workingFolder"
    return 1
fi

# Setting environment to the CMSSW_6_X release where combine tool is installed
cd $cmssw_combine
cmsenv

# Returning to the original location
cd $location

# Running the fit for root file with input templates
for filePath in $(ls -1 --color=never "$workingFolder"/*/*/*.root); do
    folderName=${filePath%/*root}
    fileName=${filePath##*/}
    fileName=${fileName%.root}
    
    echo
    echo "################################################################################"
    echo "### Performing the fit for case: $fileName"
    echo "################################################################################"
    cd $folderName
    if [ ! -d $fileName ]; then
	mkdir $fileName
    fi
    # Performing the fit
    combine -M MaxLikelihoodFit --plots --saveNormalizations --saveShapes --saveWithUncertainties --saveNLL --robustFit=1 --maxFailedSteps=15 --rMin=-100 --rMax=100 -m 125 $fileName.txt --out $fileName -v1 -n test > $fileName.log
    echo "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
    echo "Fitting done with log file:"
    echo "${filePath%.root}.log"
    echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    echo
done

# Returning to the original location and release
cd $cmssw_original
cmsenv
cd $location


