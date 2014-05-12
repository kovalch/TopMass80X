if [ "$#" -ne 1 -a "$#" -ne 2 ]; 
then
    echo " "
    echo "Use the following syntax:"
    echo "hadd_on_dcache.sh <INPUT_FOLDER_ON_DCACHE> <END_OF_FILE_NAME>"
    echo "<END_OF_FILE_NAME> is optional and may hold what should be before .root in the file name"
    echo " "
else
    INPUTFOLDER=$1
    #LFNTOPFN=`edmFileUtil -d /store`
    LFNTOPFN=srm://dcache-se-cms.desy.de:8443
    DCFOLDER=`echo ${LFNTOPFN}|grep -o /pnfs.*`
    for i in $(dcls ${DCFOLDER}${INPUTFOLDER}|grep $2.root)
    do
        echo "Deleting ${LFNTOPFN}${INPUTFOLDER}/$i"
        srm-advisory-delete ${LFNTOPFN}${INPUTFOLDER}/$i
    done
fi
