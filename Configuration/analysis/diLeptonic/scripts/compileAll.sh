#!/bin/sh


TopAnalysisDir=${CMSSW_BASE}/src/TopAnalysis
analysisDir=${TopAnalysisDir}/Configuration/analysis


ztopDir=${TopAnalysisDir}/ZTopUtils
commonSourceDir=${analysisDir}/common
mainSourceDir=${analysisDir}/diLeptonic


if [ $# -ge 2 ] ; then
    echo "Usage (default installation path): $0"
    echo "Usage (installation in path <folder>): $0 <folder>"
    exit 1
elif [ $# == 1 ]; then
    installDir=$1
    if [ ! -d "${installDir}" ]; then
        echo "Specified install directory not existing: ${installDir}"
	echo "Please create directory before installation"
	exit 2
    fi
    cd ${installDir}
    commonInstallDir=${PWD}
    mainInstallDir=${PWD}
    cd -
elif [ $# == 0 ]; then
    commonInstallDir=${commonSourceDir}
    mainInstallDir=${mainSourceDir}
fi


echo
echo
echo "Compiling project ZTopUtils using scram"
echo "Source code in folder: ${ztopDir}"
echo "Installation in folder: as defined by scram"
echo "Warning about 'Invalid tool lhapdffull' can be ignored"
echo
cd ${ztopDir}
scram b -r -j8
if [ $? -eq 0 ] ; then
    echo
    echo "Compilation successful"
else
    echo
    echo "Compilation NOT successful, stopping..."
    exit 3
fi
cd -
echo
echo


echo
echo
echo "Compiling project ntupleCommon using cmake/make"
echo "Source code in folder: ${commonSourceDir}"
echo "Installation in folder: ${commonInstallDir}"
echo
if [ ! -d "${commonInstallDir}/build_ntupleCommon" ]; then
    mkdir ${commonInstallDir}/build_ntupleCommon
fi
cd ${commonInstallDir}/build_ntupleCommon
cmake -D CMAKE_INSTALL_PREFIX=${commonInstallDir}/install ${commonSourceDir}
make -j8 install
if [ $? -eq 0 ] ; then
    echo
    echo "Compilation successful"
else
    echo
    echo "Compilation NOT successful, stopping..."
    exit 4
fi
cd -
echo
echo


echo
echo
echo "Compiling project diLeptonic using cmake/make"
echo "Source code in folder: ${mainSourceDir}"
echo "Installation in folder: ${mainInstallDir}"
echo
if [ ! -d "${mainInstallDir}/build_diLeptonic" ]; then
    mkdir ${mainInstallDir}/build_diLeptonic
fi
cd ${mainInstallDir}/build_diLeptonic
cmake -D CMAKE_INSTALL_PREFIX=${mainInstallDir}/install ${mainSourceDir}
make -j8 install
if [ $? -eq 0 ] ; then
    echo
    echo "Compilation successful"
else
    echo
    echo "Compilation NOT successful, stopping..."
    exit 5
fi
cd -
echo
echo






