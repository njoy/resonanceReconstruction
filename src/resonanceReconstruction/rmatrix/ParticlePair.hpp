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
  ParticlePairID id_;
  bool incident_;

  AtomicMass reduced_;
  double massratio_;

  //! @todo store reduced mass, wave number factor and permeability factor?

  /* auxiliary functions */
  #include "resonanceReconstruction/rmatrix/ParticlePair/src/makeID.hpp"

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
   *  @brief Return the Q value associated to the transition of the incident
   *         particle pair to this particle pair
   */
  const QValue& Q() const { return this->qvalue_; }

  /**
   *  @brief Set the Q value associated to the transition of the incident
   *         particle pair to this particle pair
   *
   *  @param[in] q   the new Q value for this particle pair
   */
  void changeQ( const QValue& q ) { this->qvalue_ = q; }

  /**
   *  @brief Return the identifier of the particle pair
   */
  const ParticlePairID& pairID() const { return this->id_; }

  /**
   *  @brief Return whether or not the particle pair is the incident particle
   *         pair
   */
  bool incident() const { return this->incident_; }

  /**
  *  @brief Toggle the incident channel state
   */
  void toggleIncident() { this->incident_ = not this->incident_; }

  /**
   *  @brief Return the mass ratio M / ( m + M ) of the particle pair
   */
  double massRatio() const { return this->massratio_; }

  #include "resonanceReconstruction/rmatrix/ParticlePair/src/reducedMass.hpp"
  #include "resonanceReconstruction/rmatrix/ParticlePair/src/mandelstam.hpp"
  #include "resonanceReconstruction/rmatrix/ParticlePair/src/waveNumber.hpp"
  #include "resonanceReconstruction/rmatrix/ParticlePair/src/sommerfeldParameter.hpp"
};
