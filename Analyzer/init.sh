module use -a /afs/desy.de/group/cms/modulefiles/
module load cmssw
cmsenv
PATH=/afs/desy.de/products/root/amd64_rhel50/5.34.00/bin:$PATH
MANPATH=/afs/desy.de/products/root/amd64_rhel50/5.34.00/man:$MANPATH
PYTHONPATH=/afs/desy.de/products/root/amd64_rhel50/5.34.00/lib:$PYTHONPATH
LD_LIBRARY_PATH=/opt/d-cache/dcap/lib64:/afs/desy.de/products/root/amd64_rhel50/5.34.00/lib:/usr/lib64/perl5/5.10.0/x86_64-linux-thread-multi/CORE:$LD_LIBRARY_PATH
ROOTSYS=/afs/desy.de/products/root/amd64_rhel50/5.34.00
source setup_env_old.sh
PATH=$PATH:`pwd`
