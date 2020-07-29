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
