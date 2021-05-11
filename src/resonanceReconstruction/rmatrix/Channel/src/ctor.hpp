private:

/**
 *  @brief Constructor
 *
 *  @param[in] channelID    the unique channel identifier
 *  @param[in] reactionID   the reaction identifier associated to the channel
 *  @param[in] incident     the current incident particle pair
 *  @param[in] pair         the particle pair associated to the channel
 *  @param[in] qValue       the Q value associated with the transition from
 *                          incident to the outgoing particle pair
 *  @param[in] numbers      the channel's quantum numbers
 *  @param[in] radii        the channel radii to be used for calculation of the
 *                          wave functions
 *  @param[in] boundary     the boundary condition (often the value of the
 *                          orbital angular momentum of the channel)
 */
Channel( const ChannelID& channelID,
         const ReactionID& reactionID,
         const ParticlePair& incident,
         const ParticlePair& pair,
         const QValue& qValue,
         const ChannelQuantumNumbers& numbers,
         const ChannelRadii& radii,
         const BoundaryCondition& boundary ) :
  id_( channelID ), reaction_( reactionID ),
  incident_( incident ), pair_( pair ), q_(qValue),
  numbers_( numbers ), radii_( radii ), boundary_( boundary ),
  spinfactor_( [&] { const auto J = numbers.totalAngularMomentum();
                     const auto ia = pair.particle().spin();
                     const auto ib = pair.residual().spin();
                     return  ( 2. * J + 1. ) / ( 2. * ia + 1. )
                                             / ( 2. * ib + 1. ); }() ) {}

public:

//! @todo pybind11 variant needs default constructor workaround
#ifdef PYBIND11
/**
 *  @brief Default constructor - only enabled for pybind11
 */
Channel() = default;
#endif

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
 *  @param[in] incident   the current incident particle pair
 *  @param[in] pair       the particle pair associated to the channel
 *  @param[in] qValue     the Q value associated with the transition from
 *                        incident to the outgoing particle pair
 *  @param[in] numbers    the channel's quantum numbers
 *  @param[in] radii      the channel radii to be used for calculation of the
 *                        wave functions
 *  @param[in] boundary   the boundary condition (often the value of the
 *                        orbital angular momentum of the channel)
 */
Channel( const ParticlePair& incident,
         const ParticlePair& pair,
         const QValue& qValue,
         const ChannelQuantumNumbers& numbers,
         const ChannelRadii& radii,
         const BoundaryCondition& boundary = 0.0 ) :
  Channel( makeChannelID( pair.pairID().symbol(), numbers ),
           ReactionID( makeReactionID( incident.pairID(), pair.pairID() ) ),
           incident, pair, qValue, numbers, radii, boundary ) {}

/**
 *  @brief Constructor
 *
 *  In a few rare cases, the automatically generated channel ID is not unique
 *  (e.g. when using multiple fission channels in which the ParticlePair for
 *  both channels is a dummy fission pair). In those cases, the user has to
 *  override the ChannelId to ensure that the ChannelId is unique over all
 *  channels.
 *
 *  This constructor allows a user to provide a replacement ParticlePairId to be
 *  specified which will only be used in the generation of the ChannelID
 *  (essentially appending the quantum numbers to the ParticlePairId).
 *
 *  @param[in] incident   the current incident particle pair
 *  @param[in] pair       the particle pair associated to the channel
 *  @param[in] id         the replacement particle pair identifier
 *  @param[in] qValue     the Q value associated with the transition from
 *                        incident to the outgoing particle pair
 *  @param[in] numbers    the channel's quantum numbers
 *  @param[in] radii      the channel radii to be used for calculation of the
 *                        wave functions
 *  @param[in] boundary   the boundary condition (often the value of the
 *                        orbital angular momentum of the channel)
 */
Channel( const ParticlePair& incident,
         const ParticlePair& pair,
         const std::string& id,
         const QValue& qValue,
         const ChannelQuantumNumbers& numbers,
         const ChannelRadii& radii,
         const BoundaryCondition& boundary = 0.0 ) :
  Channel( makeChannelID( id, numbers ),
           ReactionID( makeReactionID( incident.pairID(), pair.pairID() ) ),
           incident, pair, qValue, numbers, radii, boundary ) {}
