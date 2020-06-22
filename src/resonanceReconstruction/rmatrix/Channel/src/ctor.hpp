private:

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
  id_( id ), pair_( pair ), numbers_( numbers ),
  radii_( radii ), boundary_( boundary ),
  spinfactor_( [&] { const auto J = numbers.totalAngularMomentum();
                     const auto ia = pair.particle().spin();
                     const auto ib = pair.residual().spin();
                     return  ( 2. * J + 1. ) / ( 2. * ia + 1. )
                                             / ( 2. * ib + 1. ); }() ) {}

public:

/**
 *  @brief Constructor
 *
 *  This constructor generates the unique channel ID based on the identifier of
 *  the ParticlePair and the quantum numbers of the channel as "id{l,s,Jpi}",
 *  for example "n,Pu239_e0{0,1/2,1+}".
 *
 *  This constructor can be used when this generated channel ID is unique among
 *  all channels defined for the various SpinGroup that make up the
 *  CompoundNucleus.
 *
 *  @param[in] pair       the particle pair associated to the channel
 *  @param[in] numbers    the channel's quantum numbers
 *  @param[in] radii      the channel radii to be used for calculation of the
 *                        wave functions
 *  @param[in] boundary   the boundary condition (often the value of the
 *                        orbital angular momentum of the channel)
 */
Channel( const ParticlePair& pair,
         const ChannelQuantumNumbers& numbers,
         const ChannelRadii& radii,
         const BoundaryCondition& boundary = 0.0 ) :
  Channel( makeID( pair.pairID(), numbers ), pair, numbers, radii, boundary ) {}

/**
 *  @brief Constructor
 *
 *  In a few rare cases, the automatically generated channel ID is not unique
 *  (e.g. when using multiple fission channels in which the ParticlePair for
 *  both channels is a dummy fission pair). In those cases, the user has to
 *  override the ChannelId to insure that the ChannelId is unique over all
 *  channels.
 *
 *  This constructor allows a user to provide a replacement ParticlePairId to be
 *  specified which will only be used in the generation of the ChannelID
 *  (essentially appending the quantum numbers to the ParticlePairId).
 *
 *  @param[in] pair       the particle pair associated to the channel
 *  @param[in] id         the replacement particle pair identifier
 *  @param[in] numbers    the channel's quantum numbers
 *  @param[in] radii      the channel radii to be used for calculation of the
 *                        wave functions
 *  @param[in] boundary   the boundary condition (often the value of the
 *                        orbital angular momentum of the channel)
 */
Channel( const ParticlePair& pair,
         const ParticlePairID& id,
         const ChannelQuantumNumbers& numbers,
         const ChannelRadii& radii,
         const BoundaryCondition& boundary = 0.0 ) :
  Channel( makeID( id, numbers ), pair, numbers, radii, boundary ) {}
