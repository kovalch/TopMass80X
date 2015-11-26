#include "TopMass/TopEventTree/interface/CutFlags.h"  
#include "TopMass/TopEventTree/interface/CutFlags_LinkDef.h"
#include <algorithm>
#include <iostream>
#include <vector>

CutFlags::CutFlags() : TObject()
{
  std::string strFilterLabels[] = { "Trigger" , 
					"EventCleaning" ,
					"GoodVertex" ,
					"SignalLepton" ,
					"LooseMuonVeto" ,
 					"ElectronVeto" ,
					"ConversionRejection" ,
					"1Jet" ,
					"2Jets" ,
					"3Jets" ,
					"4Jets" ,
					"BTags" ,
				};
  filterLabels_.assign( strFilterLabels , strFilterLabels + 12 );

  init();
}

CutFlags::~CutFlags(){}

//public fkts

  void CutFlags::init()
  {
	 flags_ = Clear;  
  }

  bool CutFlags::hasPassedAll() const{
	 return flags_==allTags;

  };

  bool CutFlags::hasPassedTil(Filter fi) const{  
	return hasPassedTil( (Short_t)fi );	
 };

  bool CutFlags::hasPassedTil(Short_t sh) const{
	bool ret=true;
	switch(sh){
	  case BTags:
		ret = ret && hasPassed(BTags);
	  case Jets4:
		ret = ret && hasPassed(Jets4);
	  case Jets3:
		ret = ret && hasPassed(Jets3);	
	  case Jets2:
		ret = ret && hasPassed(Jets2);	
	  case Jet1:
		ret = ret && hasPassed(Jet1);
	  case ConversionRejection:
		ret = ret && hasPassed(ConversionRejection);	
	  case ElectronVeto:
		ret = ret && hasPassed(ElectronVeto);	
	  case LooseMuonVeto:
		ret = ret && hasPassed(LooseMuonVeto);
	  case SignalLepton:
		ret = ret && hasPassed(SignalLepton);
	  case GoodVertex:
		ret = ret && hasPassed(GoodVertex);
	  case EventCleaning:
		ret = ret && hasPassed(EventCleaning);
	  case Trigger:
		ret = ret && hasPassed(Trigger);
		break;
	default:
		ret=false;		  
	}
	return ret;
  };

  bool CutFlags::hasPassedTil(std::string str){
	return hasPassedTil( 1 << searchForLabel(str) );
  };

  bool CutFlags::hasPassed(Filter fi) const{
	return  hasPassed( (Short_t)fi );
  };

  bool CutFlags::hasPassed(Short_t sh) const{
	return  (bool)(flags_ & sh);
  };

  bool CutFlags::hasPassed(std::string str){
	return  hasPassed( 1 << searchForLabel(str) );
 };

  void CutFlags::setFlag(Filter fi){
	setFlag( (Short_t)fi );
  }; 

  void CutFlags::setFlag(Short_t sh){
	flags_ |= sh;
  };

  void CutFlags::setFlag(std::string str){
	setFlag( (Short_t)( 1 << searchForLabel(str) ) );
 };

  Short_t CutFlags::getMaskFromLabel(std::string str){
	return  (Short_t)( 1 << searchForLabel(str) );
  };

  std::string CutFlags::getLabelFromMask(Short_t sh) const{
	switch(sh){

	  case BTags:
		return filterLabels_.at(11);
	  case Jets4:
		return filterLabels_.at(10);
	  case Jets3:
		return filterLabels_.at(9);
	  case Jets2:
		return filterLabels_.at(8);
	  case Jet1:
		return filterLabels_.at(7);
	  case ConversionRejection:
		return filterLabels_.at(6);
	  case ElectronVeto:
		return filterLabels_.at(5);
	  case LooseMuonVeto:
		return filterLabels_.at(4);
	  case SignalLepton:
		return filterLabels_.at(3);
	  case GoodVertex:
		return filterLabels_.at(2);
	  case EventCleaning:
		return filterLabels_.at(1);
	  case Trigger:
		return filterLabels_.at(0);
	  case 0:
		return "Clear";
	  case allTags:
		 return "allTags";
	default:
		std::cout<<"Unvalid Call of std::string CutFlags::GetLabelFromMask(Short_t s), s is not a valid Mask"<<std::endl;
		return "";		  
	}
	
	
  };
  
  std::string CutFlags::getHasPassedTil() const{
	unsigned int ij = 0 ;
	while( ij < filterLabels_.size() && hasPassed( 1 << ij) ){
			 ++ij;	
	}
	if(ij==0){ 
		std::cout<<"CutFlags::getHasPassedTil(): Trigger Filter failed"<<std::endl; 
		return "Non";
 	}
	return filterLabels_.at( ij - 1 );
  };  


  unsigned int CutFlags::getNumberOfFilterLabels() const{
	return filterLabels_.size();
  };

  std::string CutFlags::getFilterLabel(unsigned int labelnr) const{
	return filterLabels_.at(labelnr);
  };


//private fkts

  unsigned int CutFlags::searchForLabel(std::string str){
	 std::vector<std::string>::iterator li = std::find( filterLabels_.begin() , filterLabels_.end() , str );
	if(li == filterLabels_.end()){
		 std::cout<<"Cutflags Filterlabel "<<str<<" were not found, try one of the following:"<<std::endl; 
		 printAllLabels();
		std::cout<<std::endl;
	}		
	return li - filterLabels_.begin();
 };

  void CutFlags::printAllLabels(){
	for(unsigned int i =0 ; i < filterLabels_.size() ; ++i){
		std::cout<<filterLabels_.at(i)<<std::endl;
	}
 }:
