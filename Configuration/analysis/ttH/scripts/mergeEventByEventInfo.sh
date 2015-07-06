#!/bin/sh

############################################ 
#                                          #
#   Write by: Christian Contreras-Campana  #
#   email: christian.contreras@desy.de     #
#   Date: 01.07.2015                       #
#                                          #
############################################

# Description: The following script merges the event-by-event information for synchronisation files based on whether if it's signal or background process

# Useage: ./mergeEventByEventInfo.sh ttbarH125inclusive Nominal emu

# Function handles which event-by-event information files to merge
function function_name () { 

    # File hearder containing column labeling
    fileHeader="run,lumi,event,is_SL,is_DL,lep1_pt,lep1_eta,lep1_phi,lep1_iso,lep1_pdgId,lep2_pt,lep2_eta,lep2_phi,lep2_iso,lep2_pdgId,jet1_pt,jet2_pt,jet3_pt,jet4_pt,jet1_CSVv2,jet2_CSVv2,jet3_CSVv2,jet4_CSVv2,MET_pt,MET_phi,n_jets,n_btags,bWeight,ttHFCategory"

    # Internal variables
    subprocess=$1
    systematic=(Nominal JES_UP JES_DOWN)
    channel=$([ "$3" == "" ] && echo '*' || echo "$3")

    # Check process type to produce merged files
    case "${subprocess}" in

	'ttbarDilepton')
	    if [[ -z "$2" ]];then

		# Merge all systematic type
		for(( i=0; i < ${#systematic[@]}; i++ ));
		do
		    dirPath=synchronisation/${systematic[i]}"/${channel}"
		    inFile=${systematic[i]}"_${channel}_${channel}_"${subprocess}"*.csv" 
		    outFile=${systematic[i]}"_"${subprocess}".csv"

		    # Remove original file headers and empty lines
		    cat $dirPath/$inFile | sed -r "s/${fileHeader}//g" | sed '/^$/d' > $outFile  
		    sed --in-place=.tmp "1s/^/${fileHeader}\n/" $outFile
		done
	    else		
		dirPath=synchronisation/$2"/${channel}/"
                inFile=$2"_*_*_"${subprocess}"*.csv"
		
		# If 3rd argument is empty merge all 3 channels else just over specified channel
		if [[ -z "$3" ]];then
                    outFile=$2"_"${subprocess}".csv"
		else
		    outFile=$2"_${channel}_${channel}_"${subprocess}".csv"
                fi

		 # Remove original file headers and empty lines  
                cat $dirPath/$inFile | sed -r "s/${fileHeader}//g" | sed '/^$/d' > $outFile 
		sed --in-place=.tmp "1s/^/${fileHeader}\n/" $outFile
	    fi
	    ;;
	'ttbarH125inclusive')
	    if [[ -z "$2" ]];then
		
		# Merge all systematic type 
                for(( i=0; i < ${#systematic[@]}; i++ ));
                do
                    dirPath=synchronisation/${systematic[i]}"/${channel}"
                    inFile=${systematic[i]}"_${channel}_${channel}_"${subprocess}"*.csv"
                    outFile=${systematic[i]}"_"${subprocess}".csv"
		    
		     # Remove original file headers and empty lines  
                    cat $dirPath/$inFile | sed -r "s/${fileHeader}//g" | sed '/^$/d' > $outFile 
		    sed --in-place=.tmp "1s/^/${fileHeader}\n/" $outFile
                done
	    else
		dirPath=synchronisation/$2"/${channel}"
                inFile=$2"_${channel}_${channel}_"${subprocess}"*.csv"

		# If 3rd argument is empty merge all 3 channels else just over specified channel 
		if [[ -z "$3" ]];then
                    outFile=$2"_"${subprocess}".csv"
                else
                    outFile=$2"_${channel}_${channel}_"${subprocess}".csv"
                fi

		 # Remove original file headers and empty lines  
		cat $dirPath/$inFile | sed -r "s/${fileHeader}//g" | sed '/^$/d' > $outFile 
		sed --in-place=.tmp "1s/^/${fileHeader}\n/" $outFile
                #alternative -i.tmp instead of --in-place=.tmp
            fi
	    ;;
	*) echo "Parameter number $1 is not processed"
	    ;;
    esac
}

# Function call is determine by command line parameter
if [[ -z "$1" && -z "$2" && -z "$3" ]];then
    function_name 'ttbarDilepton'
    function_name 'ttbarH125inclusive'
else
    function_name $1 $2 $3
fi    

# Remove temparary files used in the merging process
rm *.csv.tmp