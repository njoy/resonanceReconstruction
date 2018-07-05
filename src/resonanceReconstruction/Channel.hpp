using ChannelID = std::string;

/**
 *  @class
 *  @brief A reaction channel
 */
class Channel {

  /* fields */
  ChannelID id_;
  const ParticlePair& pair_; // pairs can be shared among channels
  QuantumNumbers numbers_;

  //! @todo the channel has radii for P, S and phi

public:

  /* constructor */
  Channel( const ChannelID& id,
           const ParticlePair& pair,
           const QuantumNumbers& numbers ) :
    id_( id ), pair_( pair ), numbers_( numbers ) {}

  const ChannelID& channelID() const { return this->id_; }
  const ParticlePair& particlePair() const { return this->pair_; }
  const QuantumNumbers& quantumNumbers() const { return this->numbers_; }

  double statisticalSpinFactor() const {

    auto J = this->quantumNumbers().totalAngularMomentum();
    auto ia = this->ParticlePair().particle().spin();
    auto ib = this->ParticlePair().residual().spin();
    return  ( 2. * J + 1. ) / ( 2. * ia + 1. ) / ( 2. * ib + 1. );
  }

  //! @todo the channel P, S, phi as a function of energy

};
