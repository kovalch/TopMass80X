#!/bin/env python



import fileinput,os,subprocess,sys



class config:
  dcache_desy   = "dcap://dcache-cms-dcap.desy.de"
  dcache_prefix = "/pnfs/desy.de/cms/tier2"
#  dcache_cms    = "/store/user/mseidel"
  dcache_cms    = "/store/user/eschliec/Skim_03_Data01"
  L_samples     = [
                   #'TT_cluster_8TeV-sherpa',
                   #'TT_lund_8TeV-sherpa'
                   'Skim_MJP12B_v1_data',
                   'Skim_MJP12C2_v1_data',
                   'Skim_MJP12C1_v1_data',
                   'Skim_MJP12D1_v1_data',
                   'Skim_MJP12D2_v1_data',
                   'Skim_MJP12D3_v1_data',
                  ]
  dataset_file  = "datasets.dbs"
  list_file     = "tmp_files.txt"



def initialise():
  modules = os.environ["MODULE_LIST"]
  print '---> $MODULE_LIST = '+ modules
  moduleMissing = False
  if not ('cmssw' in modules):   
    print '---> cmssw module not loaded'
    moduleMissing = True
  if not ('glite' in modules):
    print '---> glite module not loaded'
    moduleMissing = True
  if not ('dctools' in modules):
    print '---> dctools module not loaded'
    moduleMissing = True
  if moduleMissing:
    sys.exit(1)

  if os.path.isfile(config.dataset_file):
    print '---> ' + config.dataset_file + ' already exists, appending'

  

def cleanUp():
   if os.path.isfile(config.list_file):
     os.remove(config.list_file)
 


#def createDatasetFile(sample):
#  ls_cmd      = 'dcls ' + config.dcache_prefix + '/' + config.dcache_cms + '/' + sample
#  L_files     = subprocess.Popen(ls_cmd, stdout=subprocess.PIPE, bufsize=0, shell=True).stdout
#
#  dataset_file = open(config.dataset_file,'a') 
#  dataset_file.writelines(("[%s]\n" % sample))
#
#  for file in L_files:
#    check_cmd = 'edmFileUtil %s/%s/%s' % (config.dcache_cms, sample, file) 
#    print check_cmd
#    output = subprocess.Popen(check_cmd, stdout=subprocess.PIPE, bufsize=0, shell=True).stdout
#    for line in output:
#      info = line.split()
#      if (len(info)==9):
#        filename = info[0]
#        nevts    = info[5]
#        dataset_file.writelines(("%s = %s\n" % (filename, nevts)))
#
#  dataset_file.writelines("\n")
#  dataset_file.close()



def createDatasetFile(sample):
  ls_cmd    = 'dcls ' + config.dcache_prefix + '/' + config.dcache_cms + '/' + sample
  L_files   = subprocess.Popen(ls_cmd, stdout=subprocess.PIPE, bufsize=0, shell=True).stdout
  list_file = open(config.list_file,'w')
  for file in L_files:
    list_file.writelines(('%s/%s/%s') % (config.dcache_cms, sample, file))
  list_file.close()

  dataset_file = open(config.dataset_file,'a') 
  dataset_file.writelines(("[%s]\n" % sample))

  check_cmd = 'edmFileUtil -F ' + config.list_file
  output = subprocess.Popen(check_cmd, stdout=subprocess.PIPE, bufsize=1, shell=True).stdout
  for line in output:
    info = line.split()
    if (len(info)==9):
      print line
      filename    = info[0]
      nevts       = info[5]
      dataset_file.writelines('%s = %s\n' % (filename, nevts))

  dataset_file.writelines("\n")
  dataset_file.close()



def main():

  initialise() 
  for sample in config.L_samples:
    cleanUp()
    createDatasetFile(sample)
  cleanUp()

  print ''
  print '--> done'
  print ''




if __name__ == "__main__":
    main()
