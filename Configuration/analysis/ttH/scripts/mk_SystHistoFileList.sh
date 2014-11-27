#!/bin/zsh



if [ $# -ge 2 ] ; then
  echo "Usage: $0"
  echo "Usage: $0 <selectionRoot directory>"
  exit 1
fi


inputDir="Plots"
outputDir="FileLists_plot_systematic"

if [[ "x$1" != "x" ]] ; then
    inputDir="$1"
fi


mkdir -p ${outputDir}
rm ${outputDir}/HistoFileList_*


for systematic in $(ls -1 --color=never ${inputDir}/); do
    foreach sample (run qcd dyee dymumu dytautau ww wz zz wtolnu single ttbarbg ttbarsignal ttbarH ttbarW ttbarZ)
        foreach channel (ee emu mumu combined)
            if [ -d ${inputDir}/$systematic/$channel ] ; then
                ls -1 ${inputDir}/$systematic/$channel/*_source.root >> ${outputDir}/HistoFileList_$systematic\_$channel.txt
            fi
        end
    end
done
