#ifndef DataFormats_FbcmDigi_SiPadAmplitudeCollectoin_H
#define DataFormats_FbcmDigi_SiPadAmplitudeCollectoin_H

///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------

#include "DataFormats/FbcmDigi/interface/SiPadDigiData.h"
#include <vector>
#include <map>
#include <utility>

class SiPadDigiDataCollection {
public:
  typedef std::vector<SiPadDigiData>::const_iterator ContainerIterator;
  typedef std::pair<ContainerIterator, ContainerIterator> Range;
  typedef std::pair<unsigned int, unsigned int> IndexRange;
  typedef std::map<unsigned int, IndexRange> Registry;
  typedef std::map<unsigned int, IndexRange>::const_iterator RegistryIterator;

  SiPadDigiDataCollection() {}

  void put(Range input, unsigned int detID);
  const Range get(unsigned int detID) const;
  const std::vector<unsigned int> detIDs() const;

private:
  std::vector<SiPadDigiData> container_;
  Registry map_;
};

#endif  
