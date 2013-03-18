TOPMASSPATH=/afs/naf.desy.de/group/cms/scratch/eschliec/Releases/Selection/CMSSW_5_3_5/src/TopMass/Analyzer
source $TOPMASSPATH/setup_env.sh
#echo "LHAPATH = " $LHAPATH
echo "LD_LIBRARY_PATH = " $LD_LIBRARY_PATH
#ldd $TOPMASSPATH/TopMass
$TOPMASSPATH/TopMass $*
