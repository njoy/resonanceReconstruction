/**
 *  @brief Constructor
 *
 *  @param[in] id         the unique channel identifier
 *  @param[in] pair       the particle pair associated to the channel
 *  @param[in] numbers    the channel's quantum numbers
 *  @param[in] radii      the channel radii to be used for calculation of the
 *                        wave functions
 *  @param[in] boundary   the boundary condition (often the value of the
 *                        orbital angular momentum of the channel)
 */
Channel( const ChannelID& id,
         const ParticlePair& pair,
         const ChannelQuantumNumbers& numbers,
         const ChannelRadii& radii,
         const BoundaryCondition& boundary ) :
  id_( id ), pair_( pair ), numbers_( numbers ), radii_( radii ),
  boundary_( boundary ) {}

