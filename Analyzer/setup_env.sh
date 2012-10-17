fs flush /afs/naf.desy.de/user/e/eschliec/wd/TopMassImproved
#df -h
if [ -z "$TMPDIR" ];
then
    export TMPDIR=`mktemp -d`
fi
echo $TMPDIR
export LHAPATH=/afs/naf.desy.de/user/e/eschliec/wd/LHAPDF/share/lhapdf/PDFsets
export LD_LIBRARY_PATH=/afs/naf.desy.de/user/e/eschliec/wd/LHAPDF/lib:$LD_LIBRARY_PATH
