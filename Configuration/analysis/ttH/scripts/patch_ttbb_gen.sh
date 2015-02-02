#!/bin/sh

workDir="$CMSSW_BASE/src/TopAnalysis/Configuration/analysis/ttH"
patchDir=$workDir/data/patch_ttbb_gen
srcDir=$workDir/src

srcFile=$srcDir/HiggsAnalysis.cc
patchFile=$patchDir/HiggsAnalysis.cc
echo "########## Patching file: $srcFile"
patch -s $@ $srcFile < $patchFile
echo

srcFile=$srcDir/load_Analysis.cc
patchFile=$patchDir/load_Analysis.cc
echo "########## Patching file: $srcFile"
patch -s $@ $srcFile < $patchFile
echo

echo "Do not forget to recompile the code"
echo "To undo the patch, run: $0 -R"
