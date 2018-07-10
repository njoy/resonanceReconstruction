struct Neutron {};
struct Photon {};
struct ChargedParticle {};
struct Fission {};

/**
 *  @typedef
 *  @brief Channel types
 *
 *  Some of the data for a given channel will actually depend on the particle
 *  type of the channel. The penetrability, shift factor, phase shift and 
 *  coulomb phase shift will actually depend on whether the channel is a
 *  neutron, photon, charged particle or fission channel. This variant will
 *  allow us to capture that distinction.
 */
using ChannelType = std::variant< Neutron, Photon, ChargedParticle, Fission >;

