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
  ParticlePair pair_;
  ChannelQuantumNumbers numbers_;
  ChannelRadii radii_;
  BoundaryCondition boundary_;

  double spinfactor_;

  /* auxiliary functions */
  #include "resonanceReconstruction/rmatrix/Channel/src/makeID.hpp"

  //! @todo the channel has radii for P, S and phi which can be shared with
  //!       other channels
  //! @todo P, S and phi can be overriden with user defined functions
  //! @todo the boundary condition may be dependent on the channel radius (see
  //!       equation d.41 in the ENDF manual)
  //! @todo the default P, S, phi functions can be overridden
  //! @todo P, S, phi only depend on rho for all channel types except charged
  //!       particle channels for which it also depends on the permeability.

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/Channel/src/ctor.hpp"

  /**
   *  @brief Return the unique channel ID
   */
  const ChannelID& channelID() const { return this->id_; }

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
   *  @brief Return the Q value for going from the incident channel to this
   *         channel
   */
  const QValue& Q() const { return this->particlePair().Q(); }

  /**
   *  @brief Return the channel boundary condition
   */
  const BoundaryCondition& boundaryCondition() const { return this->boundary_; }

  #include "resonanceReconstruction/rmatrix/Channel/src/statisticalSpinFactor.hpp"
  #include "resonanceReconstruction/rmatrix/Channel/src/penetrability.hpp"
  #include "resonanceReconstruction/rmatrix/Channel/src/shiftFactor.hpp"
  #include "resonanceReconstruction/rmatrix/Channel/src/phaseShift.hpp"
  #include "resonanceReconstruction/rmatrix/Channel/src/coulombPhaseShift.hpp"
};
