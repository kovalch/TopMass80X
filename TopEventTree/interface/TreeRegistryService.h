#ifndef TopMass_TopEventTree_TreeRegistryService_h
#define TopMass_TopEventTree_TreeRegistryService_h

/*
 * TreeRegistryService.h
 *
 *  Created on: Feb 4, 2013
 *      Author: eschliec
 */

#include <string>

#include "TTree.h"

//#include "DataFormats/Common/interface/HLTPathStatus.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/ActivityRegistry.h"

namespace edm {
  //class HLTPathStatus;
}

class TreeRegistryService {
public:
  TreeRegistryService(const edm::ParameterSet&, edm::ActivityRegistry&);
  ~TreeRegistryService();

  TBranch* Branch(const char* name, void* address, const char* leaflist, Int_t bufsize = 32000) {
    return getTree()->Branch(name, address, leaflist, bufsize);
  }
  TBranch* Branch(const char* name, char* address, const char* leaflist, Int_t bufsize = 32000) {
    return getTree()->Branch(name, address, leaflist, bufsize);
  }
  TBranch* Branch(const char* name, Long_t address, const char* leaflist, Int_t bufsize = 32000) {
    return getTree()->Branch(name, address, leaflist, bufsize);
  }
  TBranch* Branch(const char* name, int address, const char* leaflist, Int_t bufsize = 32000){
    return getTree()->Branch(name, address, leaflist, bufsize);
  }
  template <class T> TBranch *Branch(const char* name, T*  obj, Int_t bufsize = 32000, Int_t splitlevel = 99) {
    return getTree()->Branch(name, obj, bufsize, splitlevel);
  }
  template <class T> TBranch* Branch(const char* name, T** obj, Int_t bufsize = 32000, Int_t splitlevel = 99) {
    return getTree()->Branch(name, obj, bufsize, splitlevel);
  }
  template <class T> TBranch* Branch(const char* name, const char* classname, T* obj, Int_t bufsize = 32000, Int_t splitlevel = 99) {
    return getTree()->Branch(name, classname, obj, bufsize, splitlevel);
  }
  template <class T> TBranch *Branch(const char* name, const char* classname, T** obj, Int_t bufsize = 32000, Int_t splitlevel = 99) {
    return getTree()->Branch(name, classname, obj, bufsize, splitlevel);
  }
  void Fill() { fillTree_ = true; }

private:

  std::string treeName_;
  std::string treeTitle_;
  TTree* tree_;
  bool fillTree_;

  TTree* getTree();
  void fillTree(const edm::StreamContext&);
  //void fillTree(const edm::Event&, const edm::EventSetup&);
  //void fillTree(std::string const&, edm::HLTPathStatus const&);
  void resetFill(const edm::StreamContext&) { fillTree_ = false; }
  //void resetFill(const edm::EventID&, const edm::Timestamp&) { fillTree_ = false; }
};

namespace edm {
   namespace service {
    // This function is needed so that there will be only on instance
    // of this service per process when "subprocesses" are being used.
    inline
    bool isProcessWideService(TreeRegistryService const*) { return true; }
  }
}

#endif /* TopMass_TopEventTree_TreeRegistryService_h */
