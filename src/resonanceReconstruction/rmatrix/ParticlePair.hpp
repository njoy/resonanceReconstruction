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
 *  The ParticlePair class gives us access to information related to the
 *  pair of particles such as the mass ratio and the reduced mass.
 */
class ParticlePair {

  /* fields */
  std::pair< Particle, Particle > pair_;
  ParticlePairID id_;

  AtomicMass reduced_;
  double massratio_;

  /* auxiliary functions */

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
   *  @brief Return the identifier of the particle pair
   */
  const ParticlePairID& pairID() const { return this->id_; }

  /**
   *  @brief Return the mass ratio of the particle pair
   *
   *  The mass ratio of the two particles is defined as follows:
   *     ratio = mb / ( ma + mb )
   *  in which ma and mb are the atomic mass values of the particles in the
   *  particle pair.
   */
  double massRatio() const { return this->massratio_; }

  /**
   *  @brief Return the reduced mass of the particle pair
   *
   *  The reduced mass mu of the two particles is defined as follows:
   *     mu = ma * mb / ( ma + mb )
   *  in which ma and mb are the atomic mass values of the particles in the
   *  particle pair.
   */
  const AtomicMass& reducedMass() const { return this->reduced_; }
};
