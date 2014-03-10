#!/bin/sh

# A relative or absolute directory can be given as argument.
# In this directory (or the current directory if non is given),
# the trigger reports and root files are added using
# sumTriggerReports2.pl and fasthadd.pl, respectively.
# In addition, the root files (including *.root.NNN) are deleted.
# The number of cores (nCores) that is used by fasthadd.pl
# can be changed below. Please, only use Integers!

nCores=8

if [ ${1} ]; then
    if [ -d ${1} ]; then
	cd ${1}
    else
	echo Directory \"$1\" does not exist! Aborting script!
	exit 1
    fi
fi

if [ ${2} ]; then
    configName=${2}
else
    configName="analyzeTopDiffXSec"
fi

numberOfSentJobs=`grep "numberOfJobs =" ./${configName}_cfg.py | sed -r 's/numberOfJobs = //'`
numberOfReadyJobs=`ls -l *.root | wc -l`

if [ $numberOfSentJobs = $numberOfReadyJobs ]; then

    if [ -e joined.txt ]; then
	
	JOINFILE=joined.txt
	
	for fileName in $(cat $JOINFILE); do
	    
	    newFileName=`echo ${fileName} | sed -r 's/root//'`
	    
	    echo "sumTriggerReports2.pl out*.txt >| ${newFileName}txt"
	    eval `echo "sumTriggerReports2.pl out*.txt >| ${newFileName}txt"`
	    
	    echo "fasthadd.pl -j ${nCores} ${fileName} *.root"
	    eval `echo "fasthadd.pl -j ${nCores} ${fileName} *.root"`
	    
	    if [ -e ${fileName} ]; then
		
		read -p "remove old files? (y/n) " RESP
		if [ "$RESP" = "y" ]; then
		    echo "rm ${configName}_cfg-*.root"
		    for rootFile in `ls ${configName}_cfg-*.root`; do
		    
			rm ${rootFile}
		    
		    done
		else
		    echo "old files NOT removed!"
		fi
	    fi

	    if ls *.root.* &> /dev/null; then
		
		echo rm `ls *.root.*`
		rm `ls *.root.*`

	    fi
	    
	done
	
    else
	
	echo "There is no joined.txt file -> aborting!"
	
    fi

elif [ $numberOfSentJobs -gt $numberOfReadyJobs ]; then

    echo "Not enough root files, probably not all jobs done yet -> aborting!"

elif [ $((numberOfSentJobs+1)) = $numberOfReadyJobs ]; then

    echo "There is one root file too many, probably already joined -> aborting!"

else

    echo "Too many root files -> aborting"

fi
