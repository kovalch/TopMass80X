#!/bin/sh

source $(dirname `readlink -f $0`)/parallelTools.sh

basedir="Closure/";
mkdir -p $basedir

set +e
rm -r Plots
rm -r PlotsPt*
set -e


for ch in ee emu mumu; do
    ls -1 selectionRoot/Nominal/${ch}/*.root | grep -v run2012 > closure_${ch}mc.txt
    rm FileLists/HistoFileList_Nominal_${ch}.txt
    echo "selectionRoot/closure/$ch/${ch}_ttbarsignalplustau_fakerun_nominal.root" > FileLists/HistoFileList_Nominal_${ch}.txt
    cat closure_${ch}mc.txt >> FileLists/HistoFileList_Nominal_${ch}.txt
    rm closure_${ch}mc.txt

    ### Make unfolding
    w
    $HISTO -t unfold -p +hyptoppt -c ${ch}
    grep -v Bin < UnfoldingResults/Nominal/${ch}/HypToppTResults.txt | awk '{print $2"  "$8"  "$10"  "$10}' >> ${basedir}closure_result_pt_${ch}.txt

    $HISTO -t unfold -p +hyptoprapidity -c ${ch}
    grep -v Bin < UnfoldingResults/Nominal/${ch}/HypTopRapidityResults.txt | awk '{print $2"  "$8"  "$10"  "$10}' >> ${basedir}closure_result_y_${ch}.txt

    ### Make control plot
    $HISTO -t cp -p +hyptoppt -c ${ch}
    $HISTO -t cp -p +hyptoprapidity -c ${ch}
done

mv Plots ${basedir}PlotsNominal

echo "Pt results Nominal"


exit

for slope in 0.002 0.004 0.006 0.008 -0.002 -0.004 -0.006 -0.008; do

   for ch in ee emu mumu; do
        ls -1 selectionRoot/Nominal/${ch}/*.root | grep -v run2012 > closure_${ch}mc.txt
        rm FileLists/HistoFileList_Nominal_${ch}.txt
        echo "selectionRoot/closure/$ch/${ch}_ttbarsignalplustau_fakerun_pttop${slope}.root" > FileLists/HistoFileList_Nominal_${ch}.txt
        cat closure_${ch}mc.txt >> FileLists/HistoFileList_Nominal_${ch}.txt
        rm closure_${ch}mc.txt

        ### Make unfolding
        w
        $HISTO -t unfold -p +hyptoppt -c ${ch}
        grep -v Bin < UnfoldingResults/Nominal/${ch}/HypToppTResults.txt | awk '{print $2"  "$8"  "$10"  "$10}' >> ${basedir}closure_result_pt_${ch}_${slope}.txt

        ### Make control plot
        $HISTO -t cp -p +hyptoppt -c ${ch}
   done

   mv Plots ${basedir}PlotsPt$slope
done

echo "Pt resultsPt"
# for s in 0.002 0.004 0.006 0.008 -0.002 -0.004 -0.006 -0.008; do echo -n "1.
# 2.
# 3.
# 4.
# 5." | paste - trueMadPt.txt trueMadPt.txt closure_result_pt_$s.txt ; done





set +e
rm -r Plots
rm -r PlotsY*
set -e

for slope in 0.080 0.160 0.240 0.320 -0.080 -0.160 -0.240 -0.320; do

    for ch in ee emu mumu; do
        ls -1 selectionRoot/Nominal/${ch}/*.root | grep -v run2012 > closure_${ch}mc.txt
        rm FileLists/HistoFileList_Nominal_${ch}.txt
        echo "selectionRoot/closure/$ch/${ch}_ttbarsignalplustau_fakerun_ytop${slope}.root" > FileLists/HistoFileList_Nominal_${ch}.txt
        cat closure_${ch}mc.txt >> FileLists/HistoFileList_Nominal_${ch}.txt
        rm closure_${ch}mc.txt

        ### Make unfolding
        w
        $HISTO -t unfold -p +hyptoprapidity -c ${ch}
        grep -v Bin < UnfoldingResults/Nominal/${ch}/HypTopRapidityResults.txt | awk '{print $2"  "$8"  "$10"  "$10}' >> ${basedir}closure_result_y_${ch}_${slope}.txt

        ### Make control plots
        $HISTO -t cp -p +hyptoprapidity -c ${ch}
    done

    mv Plots ${basedir}PlotsY$slope
done

echo
echo
echo " Now you can run 'root -l -b -q macros/compareClosure.C++' to plot all the results in the same plot"


# 
# 
# #paste closure_result_y_-0.008.txt closure_result_y_-0.006.txt closure_result_y_-0.004.txt closure_result_y_-0.002.txt closure_result_y_0.002.txt closure_result_y_0.004.txt closure_result_y_0.006.txt closure_result_y_0.008.txt  >| closuretesty.csv
# echo "Y resultsY"
# for s in 0.080 0.160 0.240 0.320 -0.080 -0.160 -0.240 -0.320; do echo -n "1.
# 2.
# 3.
# 4.
# 5.
# 6.
# 7.
# 8." | paste - trueMadY.txt trueMadY.txt closure_result_y_$s.txt ; done




################################################################################################################################
####
###
###        Before running this part of the scripts you have to run:
###             build/load_Analysis -f tau.root -s closure --closure nominal -c ${ch}
###           where ch = ee, emu, mumu
###
###
################################################################################################################################


# for slope in nominal; do
# 
#    for ch in ee emu mumu; do
#         ls -1 selectionRoot/Nominal/${ch}/*.root | grep -v run2012 > closure_${ch}mc.txt
#         rm FileLists/HistoFileList_Nominal_${ch}.txt
#         echo "selectionRoot/closure/$ch/${ch}_ttbarsignalplustau_fakerun_${slope}.root" > FileLists/HistoFileList_Nominal_${ch}.txt
#         cat closure_${ch}mc.txt >> FileLists/HistoFileList_Nominal_${ch}.txt
#         rm closure_${ch}mc.txt
# 
#         ### Make unfolding
#         w
#         $HISTO -t unfold -p +hyptoppt -c ${ch}
#         grep -v Bin < UnfoldingResults/Nominal/${ch}/HypToppTResults.txt | awk '{print $2"  "$8"  "$10"  "$10}' >> ${basedir}closure_result_${ch}_${slope}.txt
# 
#         ### Make control plot
#         $HISTO -t cp -p +hyptoppt -c ${ch}
#    done
# 
#    mv Plots ${basedir}Plots$slope
# done