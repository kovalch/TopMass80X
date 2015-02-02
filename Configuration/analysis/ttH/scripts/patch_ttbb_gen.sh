#!/bin/sh

scriptsDir=$(dirname $0)
patchDir=$scriptsDir/patch_ttbb_gen
srcDir=$(dirname $scriptsDir)

srcFile=$srcDir/src/HiggsAnalysis.cc
patchFile=$patchDir/HiggsAnalysis.cc
echo "########## Patching file: $srcFile"
patch $@ $srcFile < $patchFile
echo

srcFile=$srcDir/src/load_Analysis.cc
patchFile=$patchDir/load_Analysis.cc
echo "########## Patching file: $srcFile"
patch $@ $srcFile < $patchFile
echo

echo "Do not forget to recompile the code"
echo "To undo the patch, run: $0 -R"
