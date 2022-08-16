private:

/**
 *  @brief Private constructor
 */
SpinGroupBase( std::vector< ReactionID >&& reactions,
               Channel< Neutron >&& incident, ResonanceTableType&& table ) :
  incident_( std::move( incident ) ), table_( std::move( table ) ),
  reactions_( std::move( reactions ) ) {}

public:

/**
 *  @brief Constructor
 *
 *  @param[in] incident   the incident channel data for this l,J pair
 *  @param[in] table      the table of resonance parameters for this l,J pair
 */
SpinGroupBase( Channel< Neutron >&& incident, ResonanceTableType&& table ) :
  SpinGroupBase( makeReactionIdentifiers( incident, table ),
                 std::move( incident ), std::move( table ) ) {}
