#ifndef NJOY_R2_RMATRIX_CHANNEL
#define NJOY_R2_RMATRIX_CHANNEL

// system includes

// other includes
#include "resonanceReconstruction/Quantity.hpp"
#include "resonanceReconstruction/rmatrix/identifiers.hpp"
#include "resonanceReconstruction/rmatrix/options.hpp"
#include "resonanceReconstruction/rmatrix/BoundaryCondition.hpp"
#include "resonanceReconstruction/rmatrix/ChannelQuantumNumbers.hpp"
#include "resonanceReconstruction/rmatrix/ParticlePair.hpp"
#include "resonanceReconstruction/rmatrix/ChannelRadii.hpp"
#include "resonanceReconstruction/rmatrix/calculatePenetrability.hpp"
#include "resonanceReconstruction/rmatrix/calculateShiftFactor.hpp"
#include "resonanceReconstruction/rmatrix/calculatePhaseShift.hpp"
#include "resonanceReconstruction/rmatrix/calculateCoulombPhaseShift.hpp"

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {

/**
 *  @class
 *  @brief A reaction channel
 *
 *  The Channel class is used to used to store or calculate all channel specific
 *  data such as the incident particle pair, the outgoing particle pair, the
 *  associated Q-value, etc.
 *
 *  The incident and outgoing particle pair identifiers determine the reaction
 *  to which the channel contributes. For example: 'n,Cl35->n,Cl35' for elastic
 *  scattering.
 *
 *  The channel should always have a unique channel identifier associated to
 *  it. In most cases, such a unique identifier can be generated from the
 *  channel's particle pair identifier and (l,s,Jpi) set of quantum numbers,
 *  for example: 'n,Cl35{0,1,1+}'.
 *
 *  Under special cases such as partial fission channels, a user defined
 *  identifier may be needed since the associated particle pair is the
 *  same (since they contribute to the same reaction), for example:
 *  'fission1{0,1,1+}' and 'fission2{0,1,1+}'
 */
template < typename ChannelType >
class Channel {

  /* fields */
  ChannelID id_;
  ReactionID reaction_;
  ParticlePair incident_;
  ParticlePair pair_;
  QValue q_;
  ChannelQuantumNumbers numbers_;
  ChannelRadii radii_;
  BoundaryCondition boundary_;

  double spinfactor_;

  /* auxiliary functions */
  #include "resonanceReconstruction/rmatrix/Channel/src/makeChannelID.hpp"
  #include "resonanceReconstruction/rmatrix/Channel/src/makeReactionID.hpp"

  //! @todo the boundary condition may be dependent on the channel radius (see
  //!       equation d.41 in the ENDF manual)

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/Channel/src/ctor.hpp"

  /**
   *  @brief Return the unique channel ID
   */
  const ChannelID& channelID() const { return this->id_; }

  /**
   *  @brief Return the identifier of the reaction to which this channel
   *         contributes
   */
  const ReactionID& reactionID() const { return this->reaction_; }

  /**
   *  @brief Return the current incident particle pair in the channel
   */
  const ParticlePair& incidentParticlePair() const { return this->incident_; }

  /**
   *  @brief Return the particle pair in the channel
   */
  const ParticlePair& particlePair() const { return this->pair_; }

  /**
   *  @brief Return the (l,s,Jpi) quantum number set of the channel
   */
  const ChannelQuantumNumbers& quantumNumbers() const { return this->numbers_; }

  /**
   *  @brief Return the channel radii
   */
  const ChannelRadii& radii() const { return this->radii_; }

  /**
   *  @brief Return the Q value for going from the current incident channel
   *         to this channel
   */
  const QValue& Q() const { return this->q_; }

  /**
   *  @brief Return the channel boundary condition
   */
  const BoundaryCondition& boundaryCondition() const { return this->boundary_; }

  /**
   *  @brief Return whether or not the channel is an incident channel
   */
  bool isIncidentChannel() const {

    return this->particlePair().pairID() ==
           this->incidentParticlePair().pairID();
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
  double statisticalSpinFactor() const { return this->spinfactor_; }

  #include "resonanceReconstruction/rmatrix/Channel/src/belowThreshold.hpp"
  #include "resonanceReconstruction/rmatrix/Channel/src/sommerfeldParameter.hpp"
  #include "resonanceReconstruction/rmatrix/Channel/src/waveNumber.hpp"
  #include "resonanceReconstruction/rmatrix/Channel/src/penetrability.hpp"
  #include "resonanceReconstruction/rmatrix/Channel/src/shiftFactor.hpp"
  #include "resonanceReconstruction/rmatrix/Channel/src/phaseShift.hpp"
  #include "resonanceReconstruction/rmatrix/Channel/src/coulombPhaseShift.hpp"
};

} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
