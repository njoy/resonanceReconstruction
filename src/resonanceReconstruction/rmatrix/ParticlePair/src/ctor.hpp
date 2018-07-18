/**
 *  @brief Constructor
 *
 *  @param[in] particle     the incident or outgoing particle
 *  @param[in] residual     the target or residual nuclide
 *  @param[in] qvalue       the reaction Q value (in eV)
 *  @param[in] reaction     a unique reaction identifier
 */
ParticlePair( const Particle& particle,
              const Particle& residual,
              const QValue& qvalue,
              const ReactionID& reaction ) :
  pair_( particle, residual ),
  qvalue_( qvalue ), reaction_( reaction ) {}

