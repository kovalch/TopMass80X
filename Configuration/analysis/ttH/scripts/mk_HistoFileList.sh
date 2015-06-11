#!/bin/zsh



if [ $# -ge 2 ] ; then
  echo "Usage: $0"
  echo "Usage: $0 mva"
  echo "Usage: $0 <selectionRoot directory>"
  exit 1
fi


inputDir="selectionRoot"
outputDir="FileLists_plot"

if [[ "x$1" != "x" ]] ; then
    inputDir="$1"
    outputDir="FileLists_$1"
fi

# FIXME: should be removed after using unique name for mvaTopJets
if [[ "$1" == mva ]] ; then
    inputDir="mvaInput"
    outputDir="FileLists_mva"
fi


mkdir -p ${outputDir}
rm ${outputDir}/HistoFileList_*


for systematic in $(ls -1 --color=never ${inputDir}/); do
    foreach channel (ee emu mumu)
        if [ -d ${inputDir}/$systematic/$channel ] ; then
            ls -1 ${inputDir}/$systematic/$channel/${channel}_*.root >> ${outputDir}/HistoFileList_$systematic\_$channel.txt
            ls -1 ${inputDir}/$systematic/$channel/${channel}_*.root >> ${outputDir}/HistoFileList_$systematic\_combined.txt
        fi
    end
done
