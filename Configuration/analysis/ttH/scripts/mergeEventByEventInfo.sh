#!/bin/sh

############################################ 
#                                          #
#   Write by: Christian Contreras-Campana  #
#   email: christian.contreras@desy.de     #
#   Date: 01.07.2015                       #
#                                          #
############################################

# Description: The following script merges the event-by-event information for synchronisation files based on whether if it's signal or background process

# Useage: ./mergeEventByEventInfo.sh -p ttbarH125inclusive -s Nominal -c emu ( -o ttH_JESUp use only if -c is not called on)


function optparse {
    USAGE="Usage: `basename $0` [-hv] [-o arg] args"
    
    while getopts hvo:s:c:p:d OPT; do
	
	case "$OPT" in
            h)
		echo $USAGE
		exit 0
		;;
            v)
		echo "`basename $0` version 0.1"
		exit 0
		;;
            o)
		OUTPUT_FILE="$OPTARG"
		;;
	    d)
		echo "parse option takes in no arguments"
		DEBUG=true
		;;
	    s)
		SYSTEMATIC="$OPTARG"
		;;
	    c)
		CHANNEL="$OPTARG"
		;;
	    p)
                PROCESS="$OPTARG"
                ;;
            \?)
                # getopts issues an error message
		echo $USAGE >&2
		echo "Invalid option: -$OPTARG" >&2
		exit 1
		;;
	esac
    done
        
    # Remove the switches we parsed above.
    shift $((OPTIND-1))
}
# EOF


# Function handles which event-by-event information files to merge
function function_name () { 
    
    # Give local variables their corresponding global variable values
    process="$PROCESS"
    systematic="$SYSTEMATIC"
    channel="$CHANNEL"

    # File hearder containing column labeling
    fileHeader="run,lumi,event,is_SL,is_DL,lep1_pt,lep1_eta,lep1_phi,lep1_iso,lep1_pdgId,lep2_pt,lep2_eta,lep2_phi,lep2_iso,lep2_pdgId,jet1_pt,jet2_pt,jet3_pt,jet4_pt,jet1_CSVv2,jet2_CSVv2,jet3_CSVv2,jet4_CSVv2,MET_pt,MET_phi,n_jets,n_btags,bWeight,ttHFCategory"
    
    # Internal variables
    systematics=(Nominal JES_UP JES_DOWN)

    if [[ "${systematic}" != "*" && "${channel}" != "*" ]];then

	# Merge based on specified systematic type
	dirPath="synchronisation/${systematic}/${channel}"
	inFile="${systematic}_${channel}_${channel}_${process}*.csv"
	outFile="${systematic}_${channel}_${channel}_${process}.csv"
	
	# Check if file(s) exist
	countNumOfFiles=$(ls $dirPath/$inFile 2> /dev/null | wc -l)

	if [[ "${countNumOfFiles}" == "0" ]];then
	    printf "WARNING: The file(s) $dirPath/$inFile do not exist."
            exit
	fi
	
        # Remove original file headers and empty lines                                                 
        cat $dirPath/$inFile | sed -r "s/${fileHeader}//g" | sed '/^$/d' > $outFile
        sed --in-place=.tmp "1s/^/${fileHeader}\n/" $outFile
	cat $outFile | awk 'NR<2{print $0;next}{print $0| "sort --field-separator=',' -k 3 -n"}' > $outFile.old
	mv $outFile.old $outFile
		
    elif [[ "${systematic}" != "*" && "${channel}" == "*" ]];then

	dirPath="synchronisation/${systematic}/${channel}"
        inFile="${systematic}_*_*_${process}*.csv"
        outFile="${systematic}_${process}.csv"
	
	# Check if file(s) exist                                                                                          
        countNumOfFiles=$(ls $dirPath/$inFile 2> /dev/null | wc -l)
	
        if [[ "${countNumOfFiles}" == "0" ]];then
            prinf "WARNING: The file(s) $dirPath/$inFile do not exist."
            exit
        fi
	
        # Remove original file headers and empty lines    
        cat $dirPath/$inFile | sed -r "s/${fileHeader}//g" | sed '/^$/d' > $outFile
        sed --in-place=.tmp "1s/^/${fileHeader}\n/" $outFile
	cat $outFile | awk 'NR<2{print $0;next}{print $0| "sort --field-separator=',' -k 3 -n"}' > $outFile.old

	if [[ ! -z $OUTPUT_FILE ]];then
	    mv $outFile.old  $OUTPUT_FILE.csv
	    rm $outFile
	else
            mv $outFile.old $outFile
        fi

    elif [[ "$systematic" == "*" && "$channel" != "*" ]];then

	# Merge all systematic type
	for systematic in "${systematics[@]}" ; do

            dirPath="synchronisation/${systematic}/${channel}"
            inFile="${systematic}_${channel}_${channel}_${process}*.csv"
            outFile="${systematic}_${channel}_${channel}_${process}.csv"	    

	    # Check if file(s) exist                                                                                         
            countNumOfFiles=$(ls $dirPath/$inFile 2> /dev/null | wc -l)

            if [[ "${countNumOfFiles}" == "0" ]];then
                printf "WARNING: The file(s) $dirPath/$inFile do not exist."
                exit
            fi

            # Remove original file headers and empty lines                                        
            cat $dirPath/$inFile | sed -r "s/${fileHeader}//g" | sed '/^$/d' > $outFile
            sed --in-place=.tmp "1s/^/${fileHeader}\n/" $outFile
	    cat $outFile | awk 'NR<2{print $0;next}{print $0| "sort --field-separator=',' -k 3 -n"}' > $outFile.old
	    mv $outFile.old $outFile
	done
	
    elif [[ "${systematic}" == "*" && "${channel}" == "*" ]];then

	# Merge all systematic type
	for systematic in "${systematics[@]}" ; do

	    dirPath="synchronisation/${systematic}/${channel}"
            inFile="${systematic}_*_*_${process}*.csv"
            outFile="${systematic}_${process}.csv"
	    
	    # Check if file(s) exist                                                                                         
            countNumOfFiles=$(ls $dirPath/$inFile 2> /dev/null | wc -l)
	    
            if [[ "${countNumOfFiles}" == "0" ]];then
                prinf "WARNING: The file(s) $dirPath/$inFile do not exist."
                exit
            fi
	    
            # Remove original file headers and empty lines  
            cat $dirPath/$inFile | sed -r "s/${fileHeader}//g" | sed '/^$/d' > $outFile 
	    sed --in-place=.tmp "1s/^/${fileHeader}\n/" $outFile
	    cat $outFile | awk 'NR<2{print $0;next}{print $0| "sort --field-separator=',' -k 3 -n"}' > $outFile.old
	    mv $outFile.old $outFile	   
        done
    fi
}

# Parse command line
optparse $@ 

# If not filled then provide wildcard 
SYSTEMATIC=$([ "${SYSTEMATIC}" == '' ] && echo "*" || echo "${SYSTEMATIC}")
CHANNEL=$([ "${CHANNEL}" == '' ] && echo "*" || echo "${CHANNEL}")

# Function call is based on process type
if [[ -z "${PROCESS}" ]];then

    PROCESS="ttbarH125inclusive"
    function_name $PROCESS $SYSTEMATIC $CHANNEL

    PROCESS="ttbarDilepton"
    function_name $PROCESS $SYSTEMATIC $CHANNEL
    
else
    function_name $PROCESS $SYSTEMATIC $CHANNEL
fi

# Remove tmp files
rm *.csv.tmp 
