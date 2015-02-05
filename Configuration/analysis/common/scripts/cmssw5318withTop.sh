#!/bin/sh

#####################################################################################
##  If while compiling you get a virtual memory error please execute the next steps:
##
##    cd ${CMSSW_BASE}/src
##    ulimit -v 10000000
##    ulimit -m 10000000
##    scram b distclean
##    scram b -j8
##
#####################################################################################





###### CMS Software Version and Scram Architecture ########
CMS_version=CMSSW_5_3_18
export SCRAM_ARCH=slc6_amd64_gcc472





###### Function installing Hamburg TOP package #####
topAnalysis () {
    echo
    if [ -z "${TOP_TAG}" ]; then
        echo "Installing the HEAD version of the TopAnalysis code from GIT."
        echo "  To use specific tag, 'export TOP_TAG=<TAG_NAME>' BEFORE DOWNLOADING AND RUNNING the install script (such that TAG is used for both)."
        echo "  To see available tags, execute: 'git tag'"
        cd ${CMSSW_BASE}/src
        git clone https://$1@git.cern.ch/reps/TopAnalysis
        if [ $? -eq 0 ]; then
            echo "Successful download from GIT"
            echo
        else
            echo "Error in git clone! stopping..."
            exit 4
        fi
        cd -
        echo
    else
        echo "Installing the TopAnalysis code from GIT with tag: ${TOP_TAG}"
        cd ${CMSSW_BASE}/src
        git clone https://$1@git.cern.ch/reps/TopAnalysis
        if [ $? -eq 0 ]; then
            echo "Successful download from GIT"
            echo
        else
            echo "Error in git clone! stopping..."
            exit 5
        fi
        cd -
        cd $CMSSW_BASE/src/TopAnalysis
        git checkout ${TOP_TAG}
        cd -
    fi
    echo
}





###### Steering parameter for minimal or full installation ######
minimalInstall="False"
if [[ $# -eq 0 ]] || [[ $# -ge 3 ]] || [[ $# == 2 && $2 != "min" ]] ; then
    echo "Usage for full installation: $0 <CERN_USERNAME>"
    echo "Usage for minimal installation (analysis on nTuple level): $0 <CERN_USERNAME> min"
    exit 1
elif [ $# == 2 ] ; then
    if [ $2 == "min" ] ; then
        minimalInstall="True"
    fi
fi





###### File where to store the GIT-cloned release parts (avoid spamming AFS) ######
if [ `hostname | grep "nafhh"` ]; then
    export CMSSW_GIT_REFERENCE=/nfs/dust/cms/user/$USER/.cmsgit-cache
else
    export CMSSW_GIT_REFERENCE=/data/group/top/cmsswGitReference/.cmsgit-cache
fi





###### Setting up the release ######
scram p -s CMSSW ${CMS_version}
if [ $? -ne 0 ]; then
    echo "Error in scram while setting up the release! stopping..."
    exit 2
fi

cd ${CMS_version}/src
echo export SCRAM_ARCH=${SCRAM_ARCH}
eval `scram runtime -sh`
currentDir=`pwd -P`
if [ "$CMSSW_BASE/src" != "$currentDir" ]; then
    echo "Error in setting environment variable CMSSW_BASE, did you try to install in a subfolder of another CMSSW_X_Y_Z? stopping..."
    exit 3
fi
cd -





###### If parameter set, running the minimal installation: Only install our TopAnalysis ######
if [[ "$minimalInstall" == True ]] ; then
    topAnalysis $1
    echo "Minimal installation successfully done"
    exit 0
fi





###### PAT ######
### From: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATReleaseNotes52X (revision 186), for CMSSW_5_3_20

cd $CMSSW_BASE/src
git cms-addpkg PhysicsTools/PatAlgos
cd -


###### Electron MVA ID ######
### From: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideTopRefEventSel (revision 68), for CMSSW_5_3_18

cd $CMSSW_BASE/src
git cms-addpkg EgammaAnalysis/ElectronTools
cd -
cd $CMSSW_BASE/src/EgammaAnalysis/ElectronTools/data/
cat download.url | xargs wget
cd -


###### TQAF ######
### From: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideTQAFRecipes (revision 230), for CMSSW_5_3_18

cd $CMSSW_BASE/src
git cms-addpkg TopQuarkAnalysis/TopEventProducers
git cms-addpkg AnalysisDataFormats/TopObjects
git cms-addpkg TopQuarkAnalysis/Configuration
git cms-addpkg TopQuarkAnalysis/Examples
git cms-addpkg TopQuarkAnalysis/TopEventSelection
git cms-addpkg TopQuarkAnalysis/TopHitFit
git cms-addpkg TopQuarkAnalysis/TopJetCombination
git cms-addpkg TopQuarkAnalysis/TopKinFitter
git cms-addpkg TopQuarkAnalysis/TopObjectResolutions
git cms-merge-topic cms-analysis-tools:5_3_18-updateTopRefSel
cd -

######### MET Phi corrections ######
cd $CMSSW_BASE/src
git cms-addpkg JetMETCorrections/Type1MET
cd -

######### MVA MET ######
cd $CMSSW_BASE/src
#git cms-merge-topic cms-analysis-tools:5_3_14-updateSelectorUtils
#git cms-merge-topic cms-analysis-tools:5_3_13_patch2-testNewTau
git cms-merge-topic -u TaiSakuma:53X-met-131120-01
git-cms-merge-topic -u cms-met:53X-MVaNoPuMET-20131217-01
cd -

###### For full memory option of LHAPDF, we NEED to compile ElectroWeakAnalysis/Utilities after scram setup lhapdffull for speeding it up.
### For more information check the ElectroWeakAnalysis/Utilities/README file
cd $CMSSW_BASE/src
scram setup lhapdffull
git cms-addpkg ElectroWeakAnalysis/Utilities
cd -





###### Install our TopAnalysis ######
topAnalysis $1





##### Fix to avoid crash on MVA Met producer when no PV exist
## The JetMet people were contacted and they will fix (hopefully) the bug soon
## in the meantime...
cp $CMSSW_BASE/src/TopAnalysis/Configuration/analysis/common/hacks/RecoMET_METPUSubtraction_plugins_PFMETProducerMVA.cc $CMSSW_BASE/src/RecoMET/METPUSubtraction/plugins/PFMETProducerMVA.cc
cp $CMSSW_BASE/src/TopAnalysis/Configuration/analysis/common/hacks/RecoMET_METAnalyzers_BuildFile.xml $CMSSW_BASE/src/RecoMET/METAnalyzers/BuildFile.xml





###### Compile everything ######
#checkdeps -a

cd $CMSSW_BASE/src
scram b -j 8
if [ $? -ne 0 ]; then
    echo "Error in scram while compiling! stopping..."
    exit 6
fi
cd -

echo
echo "Full installation successfully done"

