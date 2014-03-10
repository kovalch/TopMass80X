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
CMS_version=CMSSW_5_3_11
export SCRAM_ARCH=slc5_amd64_gcc462





###### Function installing Hamburg TOP package #####
topAnalysis () {
    echo
    if [ -z "${TOP_TAG}" ]; then
        echo "Installing the HEAD version of the TopAnalysis code from GIT."
        echo "  To use specific tag, 'export TOP_TAG=<TAG_NAME>' BEFORE DOWNLOADING AND RUNNING the install script (such that TAG is used for both)."
        echo "  To see available tags, execute: 'git tag'"
        cd ${CMS_version}/src
        git clone https://git.cern.ch/reps/TopAnalysis
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
        cd ${CMS_version}/src
        git clone https://git.cern.ch/reps/TopAnalysis
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
if [ $# -ge 1 ] ; then
    if [ $# -ge 2 ] || [ $1 != "min" ] ; then
        echo "Usage for full installation: $0"
        echo "Usage for minimal installation (analysis on nTuple level): $0 min"
        exit 1
    else
        minimalInstall="True"
    fi
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
    topAnalysis
    echo "Minimal installation successfully done"
    exit 0
fi





###### PAT #####
### From: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATReleaseNotes52X#V08_09_59_CMSSW_5_3_11   (revision 148)

cd $CMSSW_BASE/src
addpkg DataFormats/PatCandidates V06-05-06-12
addpkg PhysicsTools/PatAlgos
cd -


###### Jet Energy Corrections #####



###### Electron ID #####

# Electron mva id stuff (following top reference twiki, and TQAF TWiki page rev.223)
# The following line is not working anymore, since the -d option returns some error, thus use workaround with the 5 lines after
#cvs co -r V00-00-30-01 -d EGamma/EGammaAnalysisTools UserCode/EGamma/EGammaAnalysisTools
cd $CMSSW_BASE/src
cvs co -r V00-00-30-01 UserCode/EGamma/EGammaAnalysisTools
mv UserCode/EGamma EGamma
rm -rf UserCode
cd -
cd $CMSSW_BASE/src/EGamma/EGammaAnalysisTools/data
cat download.url | xargs wget
cd -


###### ParticleFlow #####

cd $CMSSW_BASE/src
cvs co -r V15-02-06 RecoParticleFlow/PFProducer
cd -


###### TQAF #####

### From: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideTQAFRecipes#CMSSW_5_3_X  (revision 224)
cd $CMSSW_BASE/src
addpkg AnalysisDataFormats/TopObjects
addpkg TopQuarkAnalysis/Configuration
addpkg TopQuarkAnalysis/Examples
addpkg TopQuarkAnalysis/TopEventProducers
addpkg TopQuarkAnalysis/TopEventSelection
addpkg TopQuarkAnalysis/TopHitFit
addpkg TopQuarkAnalysis/TopJetCombination
addpkg TopQuarkAnalysis/TopKinFitter
addpkg TopQuarkAnalysis/TopObjectResolutions
addpkg TopQuarkAnalysis/TopSkimming V07-01-04
addpkg TopQuarkAnalysis/TopTools V06-07-13
cd -


# for full memory option of LHAPDF, we NEED to compile ElectroWeakAnalysis/Utilities after scram setup lhapdffull for speeding it up.
#For more information check the ElectroWeakAnalysis/Utilities/README file
cd $CMSSW_BASE/src
scram setup lhapdffull
addpkg ElectroWeakAnalysis/Utilities
cd -





###### Install our TopAnalysis ######
topAnalysis





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

