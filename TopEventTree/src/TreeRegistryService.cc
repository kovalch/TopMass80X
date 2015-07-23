/*
 * TreeRegistryService.cc
 *
 *  Created on: Feb 4, 2013
 *      Author: eschliec
 */

#include "TopMass/TopEventTree/interface/TreeRegistryService.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

TreeRegistryService::TreeRegistryService(const edm::ParameterSet& cfg, edm::ActivityRegistry& r):
treeName_ (cfg.getParameter<std::string>("treeName")),
treeTitle_(cfg.getParameter<std::string>("treeTitle")),
tree_(0),
fillTree_(false)
{
  r.watchPreEvent (this, &TreeRegistryService::resetFill);
  r.watchPostEvent(this, &TreeRegistryService::fillTree);
  //r.watchPostProcessPath(this, &TreeRegistryService::fillTree);
}

TreeRegistryService::~TreeRegistryService()
{
}

TTree*
TreeRegistryService::getTree()
{
  if(!tree_){
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // book tree if it is not existing, yet
    //////////////////////////////////////////////////////////////////////////////////////////////////
    edm::Service<TFileService> fs;
    if( !fs ) throw edm::Exception( edm::errors::Configuration, "TFileService is not registered in cfg file!" );
    tree_ = fs->make<TTree>(treeName_.c_str(), treeTitle_.c_str());
  }
  return tree_;
}

void
TreeRegistryService::fillTree(const edm::StreamContext&)
{
  if(fillTree_) getTree()->Fill();
}

//void
//TreeRegistryService::fillTree(std::string const& path, edm::HLTPathStatus const& status)
//{
//  std::cout << "path ended: " << path << " " << status.accept() << std::endl;
//  getTree()->Fill();
//}

// define this as a service
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h" 
DEFINE_FWK_SERVICE( TreeRegistryService );
