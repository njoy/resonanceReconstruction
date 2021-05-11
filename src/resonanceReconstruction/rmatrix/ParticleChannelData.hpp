#ifndef NJOY_R2_RMATRIX_PARTICLECHANNELDATA
#define NJOY_R2_RMATRIX_PARTICLECHANNELDATA

// system includes
#include <variant>
#include <vector>

// other includes
#include "Log.hpp"
#include "resonanceReconstruction/Quantity.hpp"
#include "resonanceReconstruction/rmatrix/ParticleChannel.hpp"

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {

/**
 *  @class
 *  @brief Data for a single particle channel
 *
 *  This class represent all data related to a single channel; the particle
 *  channel instance, the resonance energies and reduced widths. It provides
 *  an alternate method of constructing the SpinGroup and CompoundSystem for a
 *  set of channels and is most useful when the data is not presented by spin
 *  group to begin with, or if the spingroups are not unique.
 */
 class ParticleChannelData {

  /* fields */
  ParticleChannel channel_;
  std::vector< Energy > energies_;
  std::vector<  ReducedWidth > widths_;
  bool eliminated_;

  /* auxiliary functions */
  #include "resonanceReconstruction/rmatrix/ParticleChannelData/src/verifySize.hpp"

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/ParticleChannelData/src/ctor.hpp"

  /**
   *  @brief Return the channel
   */
  const ParticleChannel& channel() const { return this->channel_; }

  /**
   *  @brief Return the unique channel ID
   */
  const ChannelID& channelID() const {

    return std::visit( [] ( const auto& channel ) -> decltype(auto)
                          { return channel.channelID(); },
                       this->channel_ ); }

  /**
   *  @brief Return the reaction ID fort he reaction to which this channel
   *         constributes
   */
  const ReactionID& reactionID() const {

    return std::visit( [] ( const auto& channel ) -> decltype(auto)
                          { return channel.reactionID(); },
                       this->channel_ ); }

  /**
   *  @brief Return the current incident particle pair in the channel
   */
  const ParticlePair& incidentParticlePair() const {

    return std::visit( [] ( const auto& channel ) -> decltype(auto)
                          { return channel.incidentParticlePair(); },
                       this->channel_ ); }

  /**
   *  @brief Return the particle pair in the channel
   */
  const ParticlePair& particlePair() const {

    return std::visit( [] ( const auto& channel ) -> decltype(auto)
                          { return channel.particlePair(); },
                       this->channel_ ); }

  /**
   *  @brief Return the l,s,J,pi quantum numbers of the channel
   */
  const ChannelQuantumNumbers& quantumNumbers() const {

    return std::visit( [] ( const auto& channel ) -> decltype(auto)
                          { return channel.quantumNumbers(); },
                       this->channel_ ); }

  /**
   *  @brief Return the channel radii
   */
  const ChannelRadii& radii() const {

    return std::visit( [] ( const auto& channel ) -> decltype(auto)
                          { return channel.radii(); },
                       this->channel_ ); }

  /**
   *  @brief Return the Q value for going from the current incident channel
   *         to this channel
   */
  const QValue& Q() const {

    return std::visit( [] ( const auto& channel ) -> decltype(auto)
                          { return channel.Q(); },
                       this->channel_ );
  }

  /**
   *  @brief Return the channel boundary condition
   */
  const BoundaryCondition& boundaryCondition() const {

    return std::visit( [] ( const auto& channel ) -> decltype(auto)
                          { return channel.boundaryCondition(); },
                       this->channel_ );
  }

  /**
   *  @brief Return whether or not the channel is an incident channel
   */
  auto isIncidentChannel() const {

    return std::visit( [] ( const auto& channel )
                          { return channel.isIncidentChannel(); },
                       this->channel_ );
  }

  /**
   *  @brief Return the statistical spin factor
   *
   *  The statistical spin factor g of a channel is defined as follows:
   *     g = ( 2 * J + 1 ) / ( 2 * ia + 1 ) / ( 2 * ib + 1 )
   *  in which J is the total angular momentum of the channel and ia and ib
   *  are the spins of the particles in the particle pair.
   *
   *  For a particle pair involving a neutron (for which the particle spin is
   *  0.5), this reduces to:
   *    g = ( 2 * J + 1 ) / ( 2 * I + 1 ) / 2
   *  in which I is the target nucleus spin value.
   */
  auto statisticalSpinFactor() const {

    return std::visit( [] ( const auto& channel )
                          { return channel.statisticalSpinFactor(); },
                       this->channel_ );
  }

  /**
   *  @brief Return whether or not the energy is below the threshold for this
   *         channel
   *
   *  The incident energy is below the threshold energy for the channel if
   *      energy * ratio + q < 0.0
   *  where energy is the incident energy, ratio is the mass ratio M / ( m + M )
   *  for the incident particle pair and q is the Q value for this channel.
   *
   *  @param[in] energy   the energy to be tested
   */
  bool belowThreshold( const Energy& energy ) const {

    return std::visit( [&] ( const auto& channel )
                           { return channel.belowThreshold( energy ); },
                       this->channel_ );
  }

  /**
   *  @brief Return the value of Sommerfeld parameter for the channel at a
   *         given energy
   *
   *  The Sommerfeld parameter eta is an energy dependent quantity defined as
   *  follows:
   *     eta = z * Z * mu / ( 4 * pi * epsilon0 * hbar^2 * k )
   *  in which z and Z are the electrical charge of the particles in the
   *  particle pair, mu is the reduced mass of the particle pair, hbar is the
   *  Planck constant, k is the wave number and epsilon0 is the vacuum
   *  permittivity.
   *
   *  It is a dimensionless parameter.
   *
   *  @param energy   the energy at which the sommerfeld parameter needs to be
   *                  evaluated
   */
  double sommerfeldParameter( const Energy& energy ) const {

    return std::visit( [&] ( const auto& channel )
                           { return channel.sommerfeldParameter( energy ); },
                       this->channel_ );
  }

  /**
   *  @brief Return the wave number of the channel at a given energy
   *
   *  The wave number k is an energy dependent quantity defined as follows:
   *     hbar^2 k^2 = 2 * mu * ( energy * ratio + q )
   *  in which mu is the reduced mass of the channel's particle pair and ratio
   *  is the mass ratio M / ( m + M ) for the incident particle pair, q is the
   *  Q value associated to the transition of the incident particle pair to the
   *  channel's particle pair and hbar is the Planck constant.
   *
   *  @param energy   the energy at which the wave number needs to be evaluated
   */
  WaveNumber waveNumber( const Energy& energy ) const {

    return std::visit( [&] ( const auto& channel )
                           { return channel.waveNumber( energy ); },
                       this->channel_ );
  }

  /**
   *  @brief Return the penetrability for this channel as a function of energy
   *
   *  @param[in] energy   the energy at which the penetrability is needed
   */
  double penetrability( const Energy& energy ) const {

    return std::visit( [&] ( const auto& channel )
                           { return channel.penetrability( energy ); },
                       this->channel_ );
  }

  /**
   *  @brief Return the shift factor for this channel as a function of energy
   *
   *  @param[in] energy   the energy at which the shift factor is needed
   */
  double shiftFactor( const Energy& energy ) const {

    return std::visit( [&] ( const auto& channel )
                           { return channel.shiftFactor( energy ); },
                       this->channel_ );
  }

  /**
   *  @brief Return the phase shift for this channel as a function of energy
   *
   *  @param[in] energy   the energy at which the phase shift is needed
   */
  double phaseShift( const Energy& energy ) const {

    return std::visit( [&] ( const auto& channel )
                           { return channel.phaseShift( energy ); },
                       this->channel_ );
  }

  /**
   *  @brief Return the coulomb phase shift for this channel as a function of
   *         energy
   *
   *  @param[in] energy   the energy at which the coulomb phase shift is needed
   */
  double coulombPhaseShift( const Energy& energy ) const {

    return std::visit( [&] ( const auto& channel )
                           { return channel.coulombPhaseShift( energy ); },
                       this->channel_ );
  }

  /**
   *  @brief Return whether or not this is an eliminated channel
   */
  auto isEliminatedChannel() const { return this->eliminated_; }

  /**
   *  @brief Return the resonance energies
   */
  auto energies() const {

    return ranges::cpp20::views::all( this->energies_ );
  }

  /**
   *  @brief Return the reduced resonance widths
   */
  auto widths() const {

    return ranges::cpp20::views::all( this->widths_ );
  }
};

} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
