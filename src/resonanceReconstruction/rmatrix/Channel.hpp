/**
 *  @class
 *  @brief A reaction channel
 *
 *  The Channel class is used to used to store or calculate all channel specific
 *  data.
 */
class Channel {

  /* fields */
  ChannelID id_;
  const ParticlePair& pair_; // pairs can be shared among channels
  ChannelQuantumNumbers numbers_;
  ChannelRadii radii_;
  BoundaryCondition boundary_;
  ChannelType type_;

  //! @todo the channel has radii for P, S and phi which can be shared with
  //!       other channels
  //! @todo P, S and phi can be overriden with user defined functions
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
           const ChannelRadii& radii,
           const BoundaryCondition& boundary,
           const ChannelType& type ) :
    id_( id ), pair_( pair ), numbers_( numbers ), radii_( radii ),
    boundary_( boundary ), type_( type ) {}

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
   *  @brief Return the channel boundary condition
   */
  const BoundaryCondition& boundaryCondition() const { return this->boundary_; }

  /**
   *  @brief Return the channel type
   */
  const ChannelType& type() const { return this->type_; }

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

  /**
   *  @brief Return the penetrability for this channel as a function of energy
   *
   *  @param[in] energy   the energy at which the penetrability is needed
   */
  double penetrability( const Energy& energy ) const {
    auto function = [&]( auto type ){
      const double ratio = this->particlePair().waveNumber( energy ) *
                           this->radii().penetrabilityRadius( energy );
      const unsigned int l = this->quantumNumbers().orbitalAngularMomentum();
      return calculatePenetrability< decltype( type ) >( l, ratio );
    };
    return std::visit( function, this->type_ );
  }

  /**
   *  @brief Return the shift factor for this channel as a function of energy
   *
   *  @param[in] energy   the energy at which the shift factor is needed
   */
  double shiftFactor( const Energy& energy ) const {
    auto function = [&]( auto type ){
      const double ratio = this->particlePair().waveNumber( energy ) *
                           this->radii().shiftFactorRadius( energy );
      const unsigned int l = this->quantumNumbers().orbitalAngularMomentum();
      return calculateShiftFactor< decltype( type ) >( l, ratio );
    };
    return std::visit( function, this->type_ );
  }

  /**
   *  @brief Return the phase shift for this channel as a function of energy
   *
   *  @param[in] energy   the energy at which the phase shift is needed
   */
  double phaseShift( const Energy& energy ) const {
    auto function = [&]( auto type ){
      const double ratio = this->particlePair().waveNumber( energy ) *
                           this->radii().phaseShiftRadius( energy );
      const unsigned int l = this->quantumNumbers().orbitalAngularMomentum();
      return calculatePhaseShift< decltype( type ) >( l, ratio );
    };
    return std::visit( function, this->type_ );
  }

  /**
   *  @brief Return the coulomb phase shift for this channel as a function of
   *         energy
   *
   *  @param[in] energy   the energy at which the penetrability is needed
   */
  double coulombPhaseShift( const Energy& energy ) const {
    auto function = [&]( auto type ){
      const double eta = this->particlePair().etaParameter( energy ).value;
      const unsigned int l = this->quantumNumbers().orbitalAngularMomentum();
      return calculateCoulombPhaseShift< decltype( type ) >( l, eta );
    };
    return std::visit( function, this->type_ );
  }
};
