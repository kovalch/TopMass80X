#source ./setup_env.sh
#LD_LIBRARY_PATH=/cvmfs/cms.cern.ch/slc5_amd64_gcc462/external/gcc/4.6.2/lib64/:$LD_LIBRARY_PATH
export PATH=/afs/desy.de/user/g/garbersc/xxl/af-cms/CMSSW_7_6_3_patch2/bin/slc6_amd64_gcc493:/afs/desy.de/user/g/garbersc/xxl/af-cms/CMSSW_7_6_3_patch2/external/slc6_amd64_gcc493/bin:/cvmfs/cms.cern.ch/slc6_amd64_gcc493/cms/cmssw-patch/CMSSW_7_6_3_patch2/bin/slc6_amd64_gcc493:/cvmfs/cms.cern.ch/slc6_amd64_gcc493/cms/cmssw-patch/CMSSW_7_6_3_patch2/external/slc6_amd64_gcc493/bin:/cvmfs/cms.cern.ch/slc6_amd64_gcc493/cms/cmssw/CMSSW_7_6_3/bin/slc6_amd64_gcc493:/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/llvm/3.6-kpegke2/bin:/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/gcc/4.9.3/bin:/cvmfs/cms.cern.ch/common:/cvmfs/cms.cern.ch/bin::$PATH
export LD_LIBRARY_PATH=/afs/desy.de/user/g/garbersc/xxl/af-cms/CMSSW_7_6_3_patch2/biglib/slc6_amd64_gcc493:/afs/desy.de/user/g/garbersc/xxl/af-cms/CMSSW_7_6_3_patch2/lib/slc6_amd64_gcc493:/afs/desy.de/user/g/garbersc/xxl/af-cms/CMSSW_7_6_3_patch2/external/slc6_amd64_gcc493/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc493/cms/cmssw-patch/CMSSW_7_6_3_patch2/biglib/slc6_amd64_gcc493:/cvmfs/cms.cern.ch/slc6_amd64_gcc493/cms/cmssw-patch/CMSSW_7_6_3_patch2/lib/slc6_amd64_gcc493:/cvmfs/cms.cern.ch/slc6_amd64_gcc493/cms/cmssw-patch/CMSSW_7_6_3_patch2/external/slc6_amd64_gcc493/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc493/cms/cmssw/CMSSW_7_6_3/biglib/slc6_amd64_gcc493:/cvmfs/cms.cern.ch/slc6_amd64_gcc493/cms/cmssw/CMSSW_7_6_3/lib/slc6_amd64_gcc493:/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/llvm/3.6-kpegke2/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/gcc/4.9.3/lib64:/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/gcc/4.9.3/lib
echo "LD_LIBRARY_PATH = " $LD_LIBRARY_PATH
ldd TopMass
echo "executing: ./TopMass $*"
./TopMass $*

