#!/bin/zsh

inputDir_default="Plots_diffXS_systematic_input"
inputDir=$inputDir_default
outputDir="FileLists_plot_systematic"

if [[ "x$1" != "x" ]] ; then
    inputDir="$1"
fi

if [ $# -ge 2 ] || ! [ -d ${inputDir} ] ; then
  echo "Usage: $0 <selectionRoot directory>"
  echo " Default: <$inputDir_default>"
  echo " Used:    <$inputDir>  <-- does it exist?"
  exit 1
fi


if [ -d ${outputDir} ] ; then
    rm -r ${outputDir}
fi
mkdir -p ${outputDir}

nLists=0
for systematic in $(ls -1 --color=never ${inputDir}/); do
     foreach channel (ee emu mumu combined)
         if [ -d ${inputDir}/$systematic/$channel ] ; then
             ls -1 ${inputDir}/$systematic/$channel/*_source.root >> ${outputDir}/HistoFileList_$systematic\_$channel.txt
             let "nLists++"
         fi
     end
done

echo "File lists created: $nLists"
