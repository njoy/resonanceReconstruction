/**
 *  @class
 *  @brief A reaction channel
 */
class Channel {

  /* fields */
  ChannelID id_;
  const ParticlePair& pair_; // pairs can be shared among channels
  ChannelQuantumNumbers numbers_;
  BoundaryCondition boundary_;

  //! @todo the channel is of a certain type (neutron, photon, fission, charged)
  //!       which determine the default P, S and phi functions to use
  //! @todo the channel has radii for P, S and phi which can be shared with
  //!       other channels
  //! @todo the boundary condition may be dependent on the channel radius (see
  //!       equation d.41 in the ENDF manual)
  //! @todo the default P, S, phi functions can be overridden
  //! @todo P, S, phi only depend on rho for all channel types except charged 
  //!       particle channels for which it also depends on the permeability.

public:

  //! @todo use only move semantics in the constructor?

  /* constructor */
  Channel( const ChannelID& id,
           const ParticlePair& pair,
           const ChannelQuantumNumbers& numbers,
           const BoundaryCondition& boundary ) :
    id_( id ), pair_( pair ), numbers_( numbers ), boundary_( boundary ) {}

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
   *  @brief Return the channel boundary condition
   */
  const BoundaryCondition& boundaryCondition() const { return this->boundary_; }

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
  double statisticalSpinFactor() const {

    auto J = this->quantumNumbers().totalAngularMomentum();
    auto ia = this->particlePair().particle().spin();
    auto ib = this->particlePair().residual().spin();
    return  ( 2. * J + 1. ) / ( 2. * ia + 1. ) / ( 2. * ib + 1. );
  }

  //! @todo the channel P, S, phi as a function of energy

};
