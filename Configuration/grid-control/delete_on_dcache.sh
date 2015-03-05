DIRS=(TT_cluster_8TeV-sherpa/ TT_lund_8TeV-sherpa/)

for DIR in ${DIRS[@]}
do
  find /pnfs/desy.de/cms/tier2/store/user/mseidel/${DIR} -type f  -exec echo "deleting {}" \; -exec lcg-del -b -l -T  srmv2 srm://dcache-se-cms.desy.de:8443/srm/managerv2?SFN={} \;
  srmrmdir srm://dcache-se-cms.desy.de:8443/srm/managerv2?SFN=/pnfs/desy.de/cms/tier2/store/user/mseidel/${DIR}
done
