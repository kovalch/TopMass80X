#!/bin/sh

merged=''
groupspace='/nfs/dust/cms/group/tophh'

if [ -n "$1" ]; then
    echo "Files are copied to $1"
    read -p "Please type in your desired extension for the file names (e.g. MadSpin):" EXT
    for FOLDER in `ls -d naf_*`; do
	
	if [[ "${FOLDER}" == *Elqcd* ]] || [[ "${FOLDER}" == *WW* ]] || [[ "${FOLDER}" == *WZ* ]] || [[ "${FOLDER}" == *ZZ* ]] || [[ "${FOLDER}" == *Top* ]]; then
	    merged='MergedFiles/'
	else
	    merged=''
	fi

	cd $FOLDER
	if [[ "${FOLDER}" == *JERDn* ]]; then
	    for File in ????DiffXSec*JERDown*.*t; do
		if [[ "${File}" == *.txt ]]; then
		    merged='TriggerReports/'
		fi
		cp -vui "$File" ${groupspace}/$1/JERDown/${merged}"${File//JERDown/${EXT}JERDown}"
	    done
	elif [[ "${FOLDER}" == *JERUp* ]]; then
	    for File in ????DiffXSec*JERUp*.*t; do
		if [[ "${File}" == *.txt ]]; then
		    merged='TriggerReports/'
		fi
		cp -vui "$File" ${groupspace}/$1/JERUp/${merged}"${File//JERUp/${EXT}JERUp}"
	    done
	elif [[ "${FOLDER}" == *JESDn* ]]; then
	    for File in ????DiffXSec*JESDown*.*t; do
		if [[ "${File}" == *.txt ]]; then
		    merged='TriggerReports/'
		fi
		cp -vui "$File" ${groupspace}/$1/JESDown/${merged}"${File//JESDown/${EXT}JESDown}"
	    done
	elif [[ "${FOLDER}" == *JESUp* ]]; then
	    for File in ????DiffXSec*JESUp*.*t; do
		if [[ "${File}" == *.txt ]]; then
		    merged='TriggerReports/'
		fi
		cp -vui "$File" ${groupspace}/$1/JESUp/${merged}"${File//JESUp/${EXT}JESUp}"
	    done
	elif [[ "${FOLDER}" == *Mass*Dn* ]]; then
	    for File in ????DiffXSec*TopMassDown*.*t; do
		if [[ "${File}" == *.txt ]]; then
		    merged='TriggerReports/'
		fi
		cp -vui "$File" ${groupspace}/$1/TopMassDown/${merged}"${File//TopMassDown/${EXT}TopMassDown}"
	    done
	elif [[ "${FOLDER}" == *Mass*Up* ]]; then
	    for File in ????DiffXSec*TopMassUp*.*t; do
		if [[ "${File}" == *.txt ]]; then
		    merged='TriggerReports/'
		fi
		cp -vui "$File" ${groupspace}/$1/TopMassUp/${merged}"${File//TopMassUp/${EXT}TopMassUp}"
	    done
	elif [[ "${FOLDER}" == *MatchDn* ]]; then
	    for File in ????DiffXSec*MatchDown*.*t; do
		if [[ "${File}" == *.txt ]]; then
		    merged='TriggerReports/'
		fi
		cp -vui "$File" ${groupspace}/$1/MatchDown/${merged}"${File//MatchDown/${EXT}MatchDown}"
	    done
	elif [[ "${FOLDER}" == *MatchUp* ]]; then
	    for File in ????DiffXSec*MatchUp*.*t; do
		if [[ "${File}" == *.txt ]]; then
		    merged='TriggerReports/'
		fi
		cp -vui "$File" ${groupspace}/$1/MatchUp/${merged}"${File//MatchUp/${EXT}MatchUp}"
	    done
	elif [[ "${FOLDER}" == *ScaleDn* ]]; then
	    for File in ????DiffXSec*ScaleDown*.*t; do
		if [[ "${File}" == *.txt ]]; then
		    merged='TriggerReports/'
		fi
		cp -vui "$File" ${groupspace}/$1/ScaleDown/${merged}"${File//ScaleDown/${EXT}ScaleDown}"
	    done
	elif [[ "${FOLDER}" == *ScaleUp* ]]; then
	    for File in ????DiffXSec*ScaleUp*.*t; do
		if [[ "${File}" == *.txt ]]; then
		    merged='TriggerReports/'
		fi
		cp -vui "$File" ${groupspace}/$1/ScaleUp/${merged}"${File//ScaleUp/${EXT}ScaleUp}"
	    done
	else
	    for File in ????DiffXSec*.*t; do
		if [[ "${File}" == *.txt ]]; then
		    merged='TriggerReports/'
		fi
		cp -vui "$File" ${groupspace}/$1/${merged}"${File//Summer/${EXT}Summer}"
	    done
	fi
	cd -

    done

else
    echo "no destination folder given!"
fi