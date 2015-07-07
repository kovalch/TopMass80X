source ./setup_env.sh
#LD_LIBRARY_PATH=/cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/gcc/4.7.2-cms/lib64:/cvmfs/cms.cern.ch/slc6_amd64_gcc472/cms/cmssw-patch/CMSSW_5_3_14_patch2/external/slc6_amd64_gcc472/lib:$LD_LIBRARY_PATH
echo "LD_LIBRARY_PATH = " $LD_LIBRARY_PATH
export X509_USER_PROXY=/afs/desy.de/user/s/stadie/k5-ca-proxy.pem
./Skimmer $*
