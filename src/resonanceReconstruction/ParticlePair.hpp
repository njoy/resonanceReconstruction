using QValue = double;
using ReactionID = std::string;

const double hbar = 1.0;

/**
 *  @class
 *  @brief Information for a particle pair
 *
 *  A ParticlePair represents the two particles involved in a entrance or exit 
 *  reaction channel (we assume that the reaction is a two-body reaction). The
 *  pair consists of a "small" incident or outgoing particle (e.g. a neutron, 
 *  photon, alpha, etc.) and a "larger" target or residual nucleus (e.g. H1,
 *  He4, U235, etc.).
 *
 *  The first particle will ultimately determine the penetrability, shift
 *  factor, phase shift, etc. used for the channel they are associated with
 *  during resonance reconstruction.
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
                const ReactionId& reaction ) :
    pair_( std::make_pair< Particle, Particle >( particle, residual ) ),
    qvalue_( qvalue ), reaction_( reaction ) {}

  const Particle& particle() const { return this->pair_.first; }
  const Particle& residual() const { return this->pair_.second; }
  const QValue& Q() const { return this->qvalue_; }
  const ReactionId& reaction() const { return this->reaction_; }

  double reducedMass() const {

    auto ma = this->particle().mass();
    auto mb = this->particle().mass();
    return ma * mb / ( ma + mb );
  }

  double waveNumber( double energy ) const {

    auto mu = this->reducedMass()
    auto ma = this->particle().mass();
    auto q = this->Q();
    return sqrt( 2. * mu * ( energy * mu / ma + q ) ) / hbar;
  }

  double permeability( double energy ) const {

    auto mu = this->reducedMass()
    auto k = this->waveNumber( energy );
    auto za = this->particle().charge();
    auto zb = this->particle().charge();
    return ( za * zb * mu ) / ( hbar * k ) ;
  }
};
