/*
	author: Christoph Garbers
	date: Oct. 26 2015

*/
#ifndef CutFlags_h
#define CutFlags_h

//#define G__DICTIONARY

#include "TObject.h"
#include <vector>
#include <string>

//Class to save the Filter Results of a specific cmssw _cfg analysis
class CutFlags : public TObject {
public:

  CutFlags();
  ~CutFlags();

//enums to Mask the Bits of the saved Filter Flags
  enum Filter{
	Trigger = 0x0001 ,
	EventCleaning = 0x0002 ,
	GoodVertex = 0x0004 ,
	SignalLepton = 0x0008 ,
	LooseMuonVeto = 0x0010 ,
	ElectronVeto = 0x0020 ,
	ConversionRejection = 0x0040 ,
	Jet1 = 0x0080 ,
	Jets2 = 0x0100 ,
	Jets3 = 0x0200 ,
	Jets4 = 0x0400 ,
	BTags = 0x0800 ,
		
	allTags = 0x0FFF,
	Clear = 0x0000
  };

//fkts

  //resets the Class
  void init();

  //checkes if all Filter where passed
  bool hasPassedAll() const;

  //checks if a Filter and all before were passed, Input is the enum, the concret Mask or the Filterlabel
  bool hasPassedTil(Filter fi) const; 
  bool hasPassedTil(Short_t sh) const;
  bool hasPassedTil(std::string str);

  //checks if a Filter were passed, Input is the enum, the concret Mask or the Filterlabel
  bool hasPassed(Filter fi) const;
  bool hasPassed(Short_t sh) const;
  bool hasPassed(std::string str);

  //sets the Flag of a Single Filter
  void setFlag(Filter fi);  
  void setFlag(Short_t sh);
  void setFlag(std::string str);


//connverts the String Label of a Filter into its bitmask
 Short_t getMaskFromLabel(std::string str);

 //gives the Label of a corresponding bitmask
 std::string getLabelFromMask(Short_t sh) const;

 //checks till which Filter all Filter passed, returns "Non" if the Trigger-Filter has not been passed
  std::string getHasPassedTil() const;  

  //Get the number of possible filter labels
  unsigned int getNumberOfFilterLabels() const;


 //get the filter label nr. labelnr.
  std::string getFilterLabel(unsigned int labelnr) const;

  

private:

  unsigned int searchForLabel(std::string str);

  void printAllLabels();

  // member data
  Short_t flags_;

 //Labels of the different expected filter
 std::vector<std::string> filterLabels_; 
  
  ClassDef(CutFlags,1);


};

#endif
