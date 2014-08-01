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
fi


if [[ "$1" == mva ]] ; then
    inputDir="mvaInput"
    outputDir="FileLists_mva"
fi


mkdir -p ${outputDir}
rm ${outputDir}/HistoFileList_*


systematics=(
    Nominal
    PU_UP PU_DOWN
    TRIG_UP TRIG_DOWN
    LEPT_UP LEPT_DOWN
    JER_UP JER_DOWN
    JES_UP JES_DOWN
    BTAG_UP BTAG_DOWN
    BTAG_LJET_UP BTAG_LJET_DOWN
    BTAGDISCR_BSTAT1_UP BTAGDISCR_BSTAT1_DOWN
    BTAGDISCR_BSTAT2_UP BTAGDISCR_BSTAT2_DOWN
    BTAGDISCR_LSTAT1_UP BTAGDISCR_LSTAT1_DOWN
    BTAGDISCR_LSTAT2_UP BTAGDISCR_LSTAT2_DOWN
    BTAGDISCR_CERR1_UP BTAGDISCR_CERR1_DOWN
    BTAGDISCR_CERR2_UP BTAGDISCR_CERR2_DOWN
    KIN_UP KIN_DOWN
    TOP_PT_UP TOP_PT_DOWN
)


foreach sample (run qcd dyee dymumu dytautau ww wz zz wtolnu single ttbarbg ttbarsignal ttbarH ttbarW ttbarZ)
    foreach channel (ee emu mumu)
        foreach systematic ("${systematics[@]}")
            if [ -d ${inputDir}/$systematic/$channel ] ; then
                ls -1 ${inputDir}/$systematic/$channel/${channel}_$sample*.root >> ${outputDir}/HistoFileList_$systematic\_$channel.txt
                ls -1 ${inputDir}/$systematic/$channel/${channel}_$sample*.root >> ${outputDir}/HistoFileList_$systematic\_combined.txt
            fi
        end
    end
end
