#!/bin/sh


if [ $# -ge 1 ] && ( [[ "$@" == *"-c"* ]] || [[ "$@" == *"-s"* ]] || [[ "$@" == *"-f"* ]] ) ; then
    echo "Options '-s', '-c', '-f' are not allowed. They are set within the script."
    echo "You used options: $@"
    exit 1
fi



source $(dirname `readlink -f $0`)/parallelTools.sh


systematics=(
    Nominal
    LUMI_UP LUMI_DOWN
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
    XSEC_TT2B_UP XSEC_TT2B_DOWN
    XSEC_TTCC_UP XSEC_TTCC_DOWN
    XSEC_TTOTHER_UP XSEC_TTOTHER_DOWN
    XSEC_TTH_UP XSEC_TTH_DOWN
    XSEC_TTZ_UP XSEC_TTZ_DOWN
    FRAC_TTHF_UP FRAC_TTHF_DOWN
    FRAC_TTOTHER_UP FRAC_TTOTHER_DOWN
#   KIN_UP KIN_DOWN
#   TOP_PT_UP TOP_PT_DOWN
    MASS_UP MASS_DOWN
    MATCH_UP MATCH_DOWN
    SCALE_UP SCALE_DOWN
    POWHEGHERWIG POWHEG MCATNLO
    PERUGIA11NoCR PERUGIA11
    
)


for systematic in "${systematics[@]}" ; do
    printf "\n\033[1;1mStart running on systematic: ${systematic}\033[1;m\n"
    w
    for channel in ee emu mumu combined; do
    	$HISTO -s $systematic -c $channel $@ &
    done
done

# PDF systematics
# FIXME: Currently hardcoded. Should be automatised if possible
printf "\n\033[1;1mStart running on systematic: PDF_0_CENTRAL\033[1;m\n"
w
for channel in ee emu mumu combined; do
    $HISTO -s PDF_0_CENTRAL -c $channel $@ &
done
for id in $(seq 1 26) ; do
    printf "\n\033[1;1mStart running on systematic: PDF_$id UP/DOWN\033[1;m\n"
    for channel in ee emu mumu combined; do
	w
    	$HISTO -s PDF_$id"_UP" -c $channel $@ &
    	$HISTO -s PDF_$id"_DOWN" -c $channel $@ &
    done
done



wait

if [ "$isNAF" = 1 ]; then
    echo "Please check your jobs with qstat -u $USER | grep HistoDiff"
else
    echo "Processing all variations of nominal samples finished!"
fi





