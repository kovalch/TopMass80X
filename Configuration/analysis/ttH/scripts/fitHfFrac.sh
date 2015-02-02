#!/bin/bash

# This script is using the official Higgs analysis combine tool
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideHiggsAnalysisCombinedLimit

# USAGE: fitHfFrac scripts/fitHF.sh <folderName>

location=$PWD
cmssw_original=$CMSSW_BASE/src/
arch_original=$SCRAM_ARCH


# Configuring the release that contains the combine tool
if [[ $HOST = *desy-cms* ]]; then
  cmssw_combine="/data/group/top/HiggsCombineTool/CMSSW_6_1_2/src/"
  arch_combine="slc5_amd64_gcc472"
elif [[ $HOST = *naf* ]]; then
  cmssw_combine="/nfs/dust/cms/group/topcmsdesy/HiggsCombineTool/CMSSW_7_1_13/src/"
  arch_combine="slc6_amd64_gcc481"
fi

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
export SCRAM_ARCH="$arch_combine"
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
    combine -M MaxLikelihoodFit --plots --saveNormalizations --saveShapes --saveWithUncertainties --saveNLL --robustFit=1 --maxFailedSteps=30 --rMin=-100 --rMax=100 -m 125 $fileName.txt --out $fileName -v2 -n test > $fileName.log
    echo "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
    echo "Fitting done with log file:"
    echo "${filePath%.root}.log"
    echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    echo
done

# Returning to the original location and release
cd $cmssw_original
export SCRAM_ARCH="$arch_original"
cmsenv
cd $location


