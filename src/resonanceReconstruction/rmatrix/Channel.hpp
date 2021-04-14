#ifndef NJOY_R2_RMATRIX_CHANNEL
#define NJOY_R2_RMATRIX_CHANNEL

// system includes

// other includes
#include "resonanceReconstruction/Quantity.hpp"
#include "resonanceReconstruction/rmatrix/options.hpp"
#include "resonanceReconstruction/rmatrix/BoundaryCondition.hpp"
#include "resonanceReconstruction/rmatrix/ChannelQuantumNumbers.hpp"
#include "resonanceReconstruction/rmatrix/ParticlePair.hpp"
#include "resonanceReconstruction/rmatrix/ChannelRadii.hpp"
#include "resonanceReconstruction/rmatrix/ChannelID.hpp"
#include "resonanceReconstruction/rmatrix/ReactionID.hpp"
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
 *  data.
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
   *  @brief Return the reaction ID fort he reaction to which this channel
   *         constributes
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
   *  @brief Return the l,s,J,pi quantum numbers of the channel
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
