#ifndef NJOY_R2_RMATRIX_IDENTIFIERS
#define NJOY_R2_RMATRIX_IDENTIFIERS

// system includes
#include <string>

// other includes
#include "elementary/ParticleID.hpp"
#include "elementary/ParticlePairID.hpp"
#include "elementary/ReactionType.hpp"
#include "elementary/ReactionID.hpp"

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {

  // identifiers
  using ChannelID = std::string;
  using ParticleID = elementary::ParticleID;
  using ParticlePairID = elementary::ParticlePairID;
  using ReactionChannelID = std::string;
  using ReactionType = elementary::ReactionType;
  using ReactionID = elementary::ReactionID;

} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
