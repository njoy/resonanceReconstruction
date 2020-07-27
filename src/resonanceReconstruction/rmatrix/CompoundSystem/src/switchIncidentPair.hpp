/**
 *  @brief Change the incident particle pair
 *
 *  This function will change the incident channels of the CompoundSystem to 
 *  those channels with the requested particle pair. This will allow a user to
 *  calculate observables of the compound system associated to different
 *  incident particle pairs.
 *
 *  @param[in] incident   the new incident particle pair
 */
void switchIncidentPair( const ParticlePair& incident ) {

  if ( incident.pairID() != this->groups_.front().incidentpair().pairID() ) {

    ranges::for_each( this->groups_,
                      [&] ( auto& group )
                          { group.switchIncidentPair( incident ); } );
  }
}

