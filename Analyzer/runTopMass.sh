source /afs/naf.desy.de/group/cms/scratch/eschliec/TopMassImproved/setup_env.sh
echo "LHAPATH = " $LHAPATH
echo "LD_LIBRARY_PATH = " $LD_LIBRARY_PATH
ldd /afs/naf.desy.de/group/cms/scratch/eschliec/TopMassImproved/TopMass
/afs/naf.desy.de/group/cms/scratch/eschliec/TopMassImproved/TopMass $*
