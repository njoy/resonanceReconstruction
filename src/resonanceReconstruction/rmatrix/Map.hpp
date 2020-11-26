#ifndef NJOY_R2_MAP
#define NJOY_R2_MAP

// system includes
#include <map>

// other includes

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {

  // map
  template < typename Key, typename Value > using Map = std::map< Key, Value >;

} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
