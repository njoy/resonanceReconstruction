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
  #include "resonanceReconstruction/rmatrix/ParticlePair/src/ctor.hpp"

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

  #include "resonanceReconstruction/rmatrix/ParticlePair/src/reducedMass.hpp"
  #include "resonanceReconstruction/rmatrix/ParticlePair/src/waveNumber.hpp"
  #include "resonanceReconstruction/rmatrix/ParticlePair/src/etaParameter.hpp"
};
