TOPMASSPATH=/afs/naf.desy.de/group/cms/scratch/eschliec/TopMass_hg_devel/Analyzer
source $TOPMASSPATH/setup_env.sh
#echo "LHAPATH = " $LHAPATH
echo "LD_LIBRARY_PATH = " $LD_LIBRARY_PATH
#ldd $TOPMASSPATH/TopMass
$TOPMASSPATH/TopMass $*
