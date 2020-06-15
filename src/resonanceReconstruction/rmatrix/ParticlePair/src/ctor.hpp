/**
 *  @brief Constructor
 *
 *  In some cases the identifier of the particle pair that is generated
 *  automatically using the particle identifiers makes no sense, as would be
 *  the case for fission. This constructor can be used to override the
 *  automatically generated identifier.
 *
 *  @param[in] particle   the incident or outgoing particle
 *  @param[in] residual   the target or residual nuclide
 *  @param[in] qvalue     the reaction Q value (in eV)
 *  @param[in] id         a unique identifier for the particle pair
 *  @param[in] incident   the particle pair is the incident particle pair or not
 */
ParticlePair( const Particle& particle,
              const Particle& residual,
              const QValue& qvalue,
              const ParticlePairID& id,
              bool incident = false ) :
  pair_( particle, residual ), qvalue_( qvalue ),
  id_( id ), incident_( incident ),
  reduced_( [&] { const auto ma = particle.mass();
                  const auto mb = residual.mass();
                  return ma * mb / ( ma + mb ); }() ),
  massratio_( [&] { const auto ma = particle.mass();
                    const auto mb = residual.mass();
                    return mb / ( ma + mb ); }() ) {}

/**
 *  @brief Constructor
 *
 *  The identifier of the particle pair is generated automatically using the
 *  particle identifiers (e.g. n,U235_e0) when calling this constructor.
 *
 *  @param[in] particle   the incident or outgoing particle
 *  @param[in] residual   the target or residual nuclide
 *  @param[in] qvalue     the reaction Q value (in eV)
 *  @param[in] incident   the particle pair is the incident particle pair or not
 */
ParticlePair( const Particle& particle,
              const Particle& residual,
              const QValue& qvalue,
              bool incident = false ) :
  ParticlePair( particle, residual, qvalue,
                makeID( particle.particleID(), residual.particleID() ),
                incident ) {}