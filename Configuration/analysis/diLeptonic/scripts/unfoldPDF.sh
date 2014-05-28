#!/bin/sh

echo "Warning!"
echo " - do not run any other unfolding code when unfoldPDF is running!"
echo ""
echo "Ctrl-C to cancel or any other key to continue"
read -n 1 -s
echo ""

runSpecificVariation() {
    variation="$1"
    echo "*********************************************************"
    echo "* starting to unfold variation $variation"
    echo "*********************************************************"
    for ch in ee emu mumu; do
        rm FileLists/HistoFileList_Nominal_$ch.txt
        cp HistoFileList_Nominal_$ch.txt FileLists/
        echo "selectionRoot/$variation/$ch/${ch}_ttbarsignalplustau.root" >> FileLists/HistoFileList_Nominal_$ch.txt
    done

    #move directories to avoid overwritting them
    mv -f Plots Plots_temp
    mv -f UnfoldingResults UnfoldingResults_temp
    mv -f preunfolded preunfolded_temp
    if [ -d "Plots_temp/$variation" ] ; then mv "Plots_temp/${variation}" Plots ; fi
    if [ -d "UnfoldingResults_temp/$variation" ] ; then mv "UnfoldingResults_temp/${variation}" UnfoldingResults ; fi
    if [ -d "preunfolded_temp/$variation" ] ; then mv "preunfolded_temp/${variation}" preunfolded ; fi
    mkdir -p Plots/combined
    mkdir -p UnfoldingResults
    mkdir -p preunfolded

    # calculate inclusive xsection
    install/bin/Histo -t cp -p XSec -s Nominal -c ee -c emu -c mumu
    wait
    install/bin/Histo -t cp -p XSec -s Nominal -c combined

    # now calculate differential distributions
    for plot in `awk '{print $1}' HistoList | grep Hyp`; do
    #for plot in HypTTBarMass; do
        install/bin/Histo -t unfold -p +$plot -s Nominal &
    done
    wait

    # move back directories
    mv -f Plots "Plots_temp/${variation}"
    mv -f UnfoldingResults "UnfoldingResults_temp/${variation}"
    mv -f preunfolded "preunfolded_temp/${variation}"

    mv -f Plots_temp Plots
    mv -f UnfoldingResults_temp UnfoldingResults
    mv -f preunfolded_temp preunfolded
}


scripts/mk_HistoFileList.sh

for i in ee emu mumu combined; do
    grep -v ttbarsignalplustau.root < FileLists/HistoFileList_Nominal_$i.txt >| HistoFileList_Nominal_$i.txt 
done

runSpecificVariation PDF_0_CENTRAL
for no in `seq 1 22`; do
#for no in 1; do
    for var in UP DOWN; do
        runSpecificVariation "PDF_${no}_${var}"
    done
done

scripts/mk_HistoFileList.sh

echo "Done"
