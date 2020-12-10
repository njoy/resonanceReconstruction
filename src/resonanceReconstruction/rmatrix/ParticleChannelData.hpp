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
 *  @brief Single particle channel data
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
  auto channelID() const {

    return std::visit( [] ( const auto& channel )
                          { return channel.channelID(); },
                       this->channel_ ); }

  /**
   *  @brief Return the reaction ID fort he reaction to which this channel
   *         constributes
   */
  auto reactionID() const {

   return std::visit( [] ( const auto& channel )
                         { return channel.reactionID(); },
                       this->channel_ ); }

  /**
   *  @brief Return the current incident particle pair in the channel
   */
  auto incidentParticlePair() const {

   return std::visit( [] ( const auto& channel )
                         { return channel.incidentParticlePair(); },
                       this->channel_ ); }

  /**
   *  @brief Return the particle pair in the channel
   */
  auto particlePair() const {

   return std::visit( [] ( const auto& channel )
                         { return channel.particlePair(); },
                       this->channel_ ); }

  /**
   *  @brief Return the l,s,J,pi quantum numbers of the channel
   */
  auto quantumNumbers() const {

   return std::visit( [] ( const auto& channel )
                         { return channel.quantumNumbers(); },
                       this->channel_ ); }

  /**
   *  @brief Return the channel radii
   */
  auto radii() const {

   return std::visit( [] ( const auto& channel )
                         { return channel.radii(); },
                       this->channel_ ); }

  /**
   *  @brief Return the Q value for going from the current incident channel
   *         to this channel
   */
  auto Q() const {

   return std::visit( [] ( const auto& channel )
                         { return channel.Q(); },
                       this->channel_ ); }

  /**
   *  @brief Return the channel boundary condition
   */
  auto boundaryCondition() const {

   return std::visit( [] ( const auto& channel )
                         { return channel.boundaryCondition(); },
                       this->channel_ ); }

  /**
   *  @brief Return whether or not the channel is an incident channel
   */
  auto isIncidentChannel() const {

    return std::visit( [] ( const auto& channel )
                          { return channel.isIncidentChannel(); },
                      this->channel_ ); }

  /**
   *  @brief Return the resonance energies
   */
  auto energies() const { return ranges::view::all( this->energies_ ); }

  /**
   *  @brief Return the resonance widths
   */
  auto widths() const { return ranges::view::all( this->widths_ ); }

  /**
   *  @brief Return whether or not this is an eliminated channel
   */
  auto isEliminatedChannel() const { return this->eliminated_; }
};

} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
