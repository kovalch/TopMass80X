#!/bin/sh


BIN="@CMAKE_INSTALL_PREFIX@/bin"


if [ ! ${CMSSW_BASE} ] ; then
    echo "Environment variable CMSSW_BASE not set, cannot continue!"
    echo "First do the cmsenv in your release"
    exit 6
fi


# Check correctness of arguments
if [ $# -le 4 ] || [[ "$1" != "-d" ]] || ! [ -d "$2" ] || [[ "$3" != "-j" ]] ; then
    echo "Usage: $0 -d <folderWithFilelists> -j 3 4 6"
    echo "always specify with -d folder where filelists are"
    echo "always specify with -j which categories to merge (at least two): e.g. for  3, 4, and 6"
    exit 1
fi
for i in ${@:4}; do
    case $i in
        ''|*[!0-9]*) echo "Usage: $0 <folderWithFilelists> -j 3 4 6"
        echo "always specify with -d folder where filelists are"
        echo "always specify with -j which categories to merge (at least two): e.g. for  3, 4, and 6"
        exit 1 ;;
    esac
done


# All systematics except PDF variations
# FIXME: Should better be run on the content of a folder rather than on the predefined list of systematics
systematics=(
    Nominal
    PU_UP PU_DOWN
    TRIG_UP TRIG_DOWN
    LEPT_UP LEPT_DOWN
    JER_UP JER_DOWN
    JES_UP JES_DOWN
    BTAGDISCR_BPURITY_UP BTAGDISCR_BPURITY_DOWN
    BTAGDISCR_LPURITY_UP BTAGDISCR_LPURITY_DOWN
    BTAGDISCR_BSTAT1_UP BTAGDISCR_BSTAT1_DOWN
    BTAGDISCR_BSTAT2_UP BTAGDISCR_BSTAT2_DOWN
    BTAGDISCR_LSTAT1_UP BTAGDISCR_LSTAT1_DOWN
    BTAGDISCR_LSTAT2_UP BTAGDISCR_LSTAT2_DOWN
    BTAGDISCR_CERR1_UP BTAGDISCR_CERR1_DOWN
    BTAGDISCR_CERR2_UP BTAGDISCR_CERR2_DOWN
#    KIN_UP KIN_DOWN
    TOP_PT_UP TOP_PT_DOWN
    MASS_UP MASS_DOWN
    MATCH_UP MATCH_DOWN
    SCALE_UP SCALE_DOWN
    POWHEGHERWIG POWHEG MCATNLO
    PERUGIA11NoCR PERUGIA11
)


# Loop over all systematics and channels and merge jet categories
for systematic in "${systematics[@]}"; do
    for channel in ee emu mumu; do
        #echo "$BIN/categoryMerger -c $channel -s $systematic ${@:1}"
        $BIN/categoryMerger -c $channel -s $systematic ${@:1} &
        # FIXME: wait function needed?
    done
done



# Loop over all PDF systematics and channels and merge jet categories
for channel in ee emu mumu; do
    # FIXME: Should be determined in a more intelligent way? 
    for pdfId in $(seq 1 26); do
        echo $BIN/categoryMerger -c $channel -s PDF_$pdfId"_UP" ${@:1} &
        echo $BIN/categoryMerger -c $channel -s PDF_$pdfId"_DOWN" ${@:1} &
        # FIXME: wait function needed?
    done
    echo $BIN/categoryMerger -c $channel -s PDF_0_CENTRAL ${@:1} &
done

