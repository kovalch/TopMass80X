if [ "$#" -ne 2 -a "$#" -ne 3 ]; 
then
    echo " "
    echo "Use the following syntax:"
    echo "hadd_on_dcache.sh <OUTPUT_FILE> <INPUT_FOLDER_ON_DCACHE> <END_OF_FILE_NAME>"
    echo "<END_OF_FILE_NAME> is optional and may hold what should be before .root in the file name"
    echo " "
else
    OUTPUTFILE=$1
    #INPUTFOLDER=/store/user/eschliec/TopMassTreeWriter_02_Data02/MJP12D1_v1_data
    INPUTFOLDER=$2
    #LFNTOPFN=`edmFileUtil -d /store`
    LFNTOPFN=dcap://dcache-cms-dcap.desy.de/
    DCFOLDER=`echo ${LFNTOPFN}|grep -o /pnfs.*`
    #dcls ${DCFOLDER}${INPUTFOLDER}|grep $3.root
    #for i in \`dcls ${DCFOLDER}${INPUTFOLDER} | grep $3.root\`; do echo -n "${LFNTOPFN}${INPUTFOLDER}/$i
    hadd ${OUTPUTFILE} `for i in \`dcls ${DCFOLDER}${INPUTFOLDER}|grep $3.root\`; do echo -n "${LFNTOPFN}${INPUTFOLDER}/$i "; done`
fi
