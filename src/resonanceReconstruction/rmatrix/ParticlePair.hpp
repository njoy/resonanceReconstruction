/**
 *  @class
 *  @brief Information for a particle pair
 *
 *  A ParticlePair represents the two particles involved in a entrance or exit 
 *  reaction channel (we assume that the reaction is a two-body reaction). The
 *  pair consists of a "small" incident or outgoing particle (e.g. a neutron, 
 *  photon, alpha, etc.) and a "larger" target or residual nucleus (e.g. H1,
 *  He4, U235, etc.).
 */
class ParticlePair {

  /* fields */
  std::pair< Particle, Particle > pair_;
  QValue qvalue_;
  ReactionID reaction_;

  //! @todo store reduced mass, wave number factor and permeability factor?

public:

  /* constructor */
  ParticlePair( const Particle& particle,
                const Particle& residual,
                const QValue& qvalue,
                const ReactionID& reaction ) :
    pair_( particle, residual ),
    qvalue_( qvalue ), reaction_( reaction ) {}

  /**
   *  @brief Return the first particle of the particle pair
   */
  const Particle& particle() const { return this->pair_.first; }

  /**
   *  @brief Return the second particle of the particle pair
   */
  const Particle& residual() const { return this->pair_.second; }

  /**
   *  @brief Return the Q value associated to the particle pair
   */
  const QValue& Q() const { return this->qvalue_; }

  /**
   *  @brief Return the ID of the reaction associated to the particle pair
   */
  const ReactionID& reaction() const { return this->reaction_; }

  /**
   *  @brief Return the reduced mass of the particle pair
   *
   *  The reduced mass mu of the two particles is defined as follows:
   *     mu = ma * mb / ( ma + mb )
   *  in which ma and mb are the atomic mass values of the particles in the
   *  particle pair.
   */
  AtomicMass reducedMass() const {

    auto ma = this->particle().mass();
    auto mb = this->residual().mass();
    return ma * mb / ( ma + mb );
  }

  /**
   *  @brief Return the wave number of the particle pair at a given energy
   *
   *  The wave number k is an energy dependent quantity defined as follows:
   *     hbar^2 k^2 = 2 * mu * ( energy * mu / ma + q )
   *  in which mu is the reduced mass of the particle pair, ma is the atomic
   *  mass of the first particle in the particle pair, q is the Q value 
   *  associated to the particle pair and hbar is the Planck constant.
   *
   *  @param energy   the energy at which the wave number needs to be evaluated
   */
  WaveNumber waveNumber( const Energy& energy ) const {

    auto mu = this->reducedMass();
    auto ma = this->particle().mass();
    auto q = this->Q();
    return sqrt( 2. * mu * ( energy * mu / ma + q ) ) / hbar;
  }

  /**
   *  @brief Return the value of eta for the particle pair at a given energy
   *
   *  The parameter eta is an energy dependent quantity defined as follows:
   *     eta = za * zb * mu / hbar / k
   *  in which za and zb are the electrical charge of the particles in the 
   *  particle pair, mu is the reduced mass of the particle pair, hbar is the 
   *  Planck constant and k is the wave number.
   *
   *  @param energy   the energy at which the permeability needs to be evaluated
   */
  EtaParameter etaParameter( const Energy& energy ) const {

    auto za = this->particle().charge();
    auto zb = this->residual().charge();
    auto mu = this->reducedMass();
    auto k = this->waveNumber( energy );
    return ( za * zb * mu ) / ( hbar * k ) ;
  }
};
