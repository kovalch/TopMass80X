/*
 * RandomSubsetCreator.h
 *
 *  Created on: Oct 23, 2012
 *      Author: eschliec
 */

#ifndef RANDOMSUBSETCREATOR_H_
#define RANDOMSUBSETCREATOR_H_

#include "TTree.h"

class RandomSubsetCreator {
public:
  RandomSubsetCreator();
  virtual ~RandomSubsetCreator();
  virtual TTree* CreateRandomSubset() = 0;
};

#endif /* RANDOMSUBSETCREATOR_H_ */
