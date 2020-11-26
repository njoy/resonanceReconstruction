#ifndef NJOY_R2_RMATRIX_PARTICLECHANNEL
#define NJOY_R2_RMATRIX_PARTICLECHANNEL

// system includes
#include <variant>

// other includes
#include "resonanceReconstruction/rmatrix/ChannelTypes.hpp"
#include "resonanceReconstruction/rmatrix/Channel.hpp"

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {

/**
 *  @typedef
 *  @brief Channel types
 *
 *  Some of the data for a given channel will actually depend on the particle
 *  type of the channel. The penetrability, shift factor, phase shift and
 *  coulomb phase shift will depend on whether the channel is a neutron, photon,
 *  charged particle or fission channel. This variant will allow us to capture
 *  that distinction.
 */
using ParticleChannel = std::variant< Channel< Neutron >,
                                      Channel< Photon >,
                                      Channel< ChargedParticle >,
                                      Channel< Fission > >;

} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
