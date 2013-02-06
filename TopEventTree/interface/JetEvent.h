/*
 * JetEvent.h
 *
 *  Created on: Feb 6, 2013
 *      Author: eschliec
 */

#ifndef JETEVENT_H_
#define JETEVENT_H_

#include "TObject.h"
#include "TLorentzVector.h"
#include <vector>


class JetEvent : public TObject {
public:

  JetEvent();
  void init();

  unsigned int nJet() { return jet.size(); }

  ClassDef(JetEvent,1);

  // JetEvent data

  std::vector<TLorentzVector> jet;
  std::vector<double> jetCharge;
  std::vector<double> jetFlavour;
  std::vector<double> jetCSV;


private:
};

#endif /* JETEVENT_H_ */
