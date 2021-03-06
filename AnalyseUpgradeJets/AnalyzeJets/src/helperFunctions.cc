#include <vector>
#include <algorithm>  // for sorting

#ifndef HELPER_FUNCTIONS
#define HELPER_FUNCTIONS

// Extract the median value of an array
//double Median( std::vector<double> aVec);



// Median calculating code used in L1TowerJetPUEstimator
template <typename t>
t Median( std::vector<t> aVec){

  // Order vector collection
  std::sort( aVec.begin(), aVec.end() );

  t median(0);
  int size = aVec.size();
  int halfSize = size/2;
  if( size == 0 ){
    median = 0;
  }
  else if( size == 1 ){
    median = aVec[0];
  }
  else if( size%2 == 0 ){
    // Even number of entries, take average of the values around center
    median = ( aVec[ halfSize - 1 ] + aVec[ halfSize ] ) * 0.5;
  }
  else{
    // Odd number of entries, halfSize is central element
    median = aVec[ halfSize ];
  }

  return median;
}





#endif
