///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------


#include "DataFormats/FbcmDigi/interface/SiPadDigiDataCollection.h"
#include <iostream>
#include <algorithm>



void SiPadDigiDataCollection::put(Range input, unsigned int detID) {
  // put in Digis of detID

  // store size of vector before put
  IndexRange inputRange;

  // put in PixelDigis from input
  bool first = true;

  // fill input in temporary vector for sorting
  std::vector<SiPadDigiData> temporary;
  SiPadDigiDataCollection::ContainerIterator sort_begin = input.first;
  SiPadDigiDataCollection::ContainerIterator sort_end = input.second;
  for (; sort_begin != sort_end; ++sort_begin) {
    temporary.push_back(*sort_begin);
  }
  std::sort(temporary.begin(), temporary.end());

  // iterators over input
  SiPadDigiDataCollection::ContainerIterator begin = temporary.begin();
  SiPadDigiDataCollection::ContainerIterator end = temporary.end();
  for (; begin != end; ++begin) {
    container_.push_back(*begin);
    if (first) {
      inputRange.first = container_.size() - 1;
      first = false;
    }
  }
  inputRange.second = container_.size() - 1;

  // fill map
  map_[detID] = inputRange;
}

const SiPadDigiDataCollection::Range SiPadDigiDataCollection::get(unsigned int detID) const {
  // get Digis of detID

  auto found = map_.find(detID);
  SiPadDigiDataCollection::IndexRange returnIndexRange{};
  if (found != map_.end()) {
    returnIndexRange = found->second;
  }

  SiPadDigiDataCollection::Range returnRange;
  returnRange.first = container_.begin() + returnIndexRange.first;
  returnRange.second = container_.begin() + returnIndexRange.second + 1;

  return returnRange;
}

const std::vector<unsigned int> SiPadDigiDataCollection::detIDs() const {
  // returns vector of detIDs in map

  SiPadDigiDataCollection::RegistryIterator begin = map_.begin();
  SiPadDigiDataCollection::RegistryIterator end = map_.end();

  std::vector<unsigned int> output;

  for (; begin != end; ++begin) {
    output.push_back(begin->first);
  }

  return output;
}
