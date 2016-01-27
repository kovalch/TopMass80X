#source ./setup_env.sh
#LD_LIBRARY_PATH=/cvmfs/cms.cern.ch/slc5_amd64_gcc462/external/gcc/4.6.2/lib64/:$LD_LIBRARY_PATH
export PATH=/afs/desy.de/user/g/garbersc/xxl/af-cms/CMSSW_7_4_12_patch4/bin/slc6_amd64_gcc491:/afs/desy.de/user/g/garbersc/xxl/af-cms/CMSSW_7_4_12_patch4/external/slc6_amd64_gcc491/bin:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw-patch/CMSSW_7_4_12_patch4/bin/slc6_amd64_gcc491:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw-patch/CMSSW_7_4_12_patch4/external/slc6_amd64_gcc491/bin:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_12/bin/slc6_amd64_gcc491:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/llvm/3.6/bin:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/gcc/4.9.1-cms/bin:/cvmfs/cms.cern.ch/common:/cvmfs/cms.cern.ch/bin:$PATH
export LD_LIBRARY_PATH=/afs/desy.de/user/g/garbersc/xxl/af-cms/CMSSW_7_4_12_patch4/biglib/slc6_amd64_gcc491:/afs/desy.de/user/g/garbersc/xxl/af-cms/CMSSW_7_4_12_patch4/lib/slc6_amd64_gcc491:/afs/desy.de/user/g/garbersc/xxl/af-cms/CMSSW_7_4_12_patch4/external/slc6_amd64_gcc491/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw-patch/CMSSW_7_4_12_patch4/biglib/slc6_amd64_gcc491:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw-patch/CMSSW_7_4_12_patch4/lib/slc6_amd64_gcc491:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw-patch/CMSSW_7_4_12_patch4/external/slc6_amd64_gcc491/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_12/biglib/slc6_amd64_gcc491:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_12/lib/slc6_amd64_gcc491:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/llvm/3.6/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/gcc/4.9.1-cms/lib64:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/gcc/4.9.1-cms/lib
echo "LD_LIBRARY_PATH = " $LD_LIBRARY_PATH
ldd TopMass
echo "executing: ./TopMass $*"
./TopMass $*

