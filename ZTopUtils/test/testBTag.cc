/*
 * testBTag.cc
 *
 *  Created on: Apr 8, 2015
 *      Author: kiesej
 */

#include "../interface/bTagBase.h"
#include "../interface/bTagSFBase.h"
#include "TRandom3.h"
#include "../interface/miscUtils.h"

#include "TFile.h"

class oldBtagWrapper : public ztop::bTagBase{
public:
	bool jetIsTagged(const float&pt, const float&eta, const int & genPartonFlavor,
			const float & tagValue, const unsigned int & seed) const{
		return bTagBase::jetIsTagged(pt,eta,genPartonFlavor,tagValue,seed);
	}


};

int main(){

	using namespace ztop;

	bTagSFBase::debug=true;
	oldBtagWrapper oldb;
	oldb.setWorkingPoint(oldBtagWrapper::csvm_wp);
	oldb.setMakeEff(true);
	oldb.setSampleName("bla");

	//	bTagSFBase::debug=true;

	bTagSFBase newbt;
	newbt.loadSF("../data/CSV_8TeV.csv",BTagEntry::OP_MEDIUM,"csv","comb","up","down");
	newbt.setMakeEff(true);

	bTagSFBase btnew2,btnew3,btnew4;

	//btnew3.loadSF("../CSV.csv",BTagEntry::OP_RESHAPING,"csv","comb","up","down");
	//btnew4.loadSF("../CSV.csv",BTagEntry::OP_RESHAPING,"csv","mujets","up","down");

	//fails due to missing entries in the file
	// btnew3.loadSF("../data/CSV_8TeV.csv",BTagEntry::OP_RESHAPING,"csv","comb","up","down");

	btnew2.loadUDSGSF("../data/CSV_8TeV.csv",BTagEntry::OP_MEDIUM,"csv","comb","up","down");
	btnew2.loadBCSF("../data/CSV_8TeV.csv",BTagEntry::OP_MEDIUM,"csv","mujets","up","down");



	std::cout << "filling effs" <<std::endl;
	TRandom3 * rand=new TRandom3();

	for(int flav=1;flav<=23;flav++){
		for(float eta=-2.4;eta<2.4;eta+=0.001){
			for(float i=20;i<800;i++){
				float disc=rand->Gaus(0.5,0.5);
				if(flav==4)
					disc=rand->Gaus(0.3,0.5);
				else if(flav!=5)
					disc=rand->Gaus(0.1,0.5);
				if(disc>1)disc=1;
				if(disc<0)disc=-1;
				oldb.fillEff(i,fabs(eta),flav,disc,1);
				newbt.fillEff(i,eta,flav,disc,1);

			}
		}
		if(flav==5)flav=23;
	}

	oldb.makeEffs();

	std::cout << "write out effs" <<std::endl;
	TFile * f=0;
	f=new TFile("testbtag.root","RECREATE");
	newbt.writeToTFile(f);
	f->Close();
	delete f;

	oldb.setMakeEff(false);

	f=new TFile("testbtag.root","OPEN");
	std::cout << "read in effs" <<std::endl;
	btnew2.readFromTFile(f);
	btnew3.readFromTFile(f);

	btnew2.setSystematics(bTagSFBase::lightdown);
	btnew2.setMakeEff(false);
	btnew3.setMakeEff(false);
	oldb.setSampleName("bla");
	//oldb.setSystematic(bTagBase::heavydown);
	std::cout << "apply SFs" <<std::endl;
	size_t same=0,diff=0,oldtag=0,newtag=0;


	BTagEntry::JetFlavor bjetflav=BTagEntry::FLAV_B;

	for(size_t n=0;n<1;n++){
		for(int flav=1;flav<=23;flav++){
			for(float eta=-2.41;eta<=2.39;eta+=0.01){
				for(float i=21;i<800;i++){
					float disc=rand->Gaus(0.5,0.5);
					if(flav==4)
						disc=rand->Gaus(0.3,0.5);
					else if(flav!=5)
						disc=rand->Gaus(0.1,0.5);
					if(disc>1)disc=1;
					if(disc<0)disc=-1;
					//disc=0.898;
					bool newsftag= btnew2.jetIsTagged(i,eta,flav,disc,3*i) ;
					bool oldsftag= oldb.jetIsTagged(i,fabs(eta),flav,disc,3*i) ;

					//std::cout << sf1 << " " << sf2<<std::endl;

					//	std::cout << btnew3.getJetDiscrShapeWeight(i,fabs(eta),flav,disc) << std::endl;;

					if((newsftag && oldsftag) || (!newsftag && !oldsftag))
						same++;
					else
						diff++;
					if(newsftag)
						newtag++;
					if(oldsftag)
						oldtag++;
				}
			}
			if(flav==5)flav=23;
		}
	}
	ZTOP_COUTVAR(same);
	ZTOP_COUTVAR(diff);
	ZTOP_COUTVAR(oldtag);
	ZTOP_COUTVAR(newtag);
	float reldifftagged=(((float)newtag-(float)oldtag)/(float)oldtag);
	ZTOP_COUTVAR(reldifftagged);
	float reldiffdec=((float)diff/(float)(diff+same));
	ZTOP_COUTVAR(reldiffdec);
	//std::cout << same << " " << diff << " " << oldtag << " " << newtag << std::endl;
	delete rand;
}


