source ./setup_env.sh
LD_LIBRARY_PATH=/cvmfs/cms.cern.ch/slc5_amd64_gcc462/external/gcc/4.6.2/lib64/:$LD_LIBRARY_PATH
echo "LD_LIBRARY_PATH = " $LD_LIBRARY_PATH
ldd TopMass
echo "executing: ./TopMass $*"
./TopMass $*

