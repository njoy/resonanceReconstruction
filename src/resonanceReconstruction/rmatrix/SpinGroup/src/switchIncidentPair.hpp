/**
 *  @brief Change the incident particle pair
 *
 *  This function will change the incident channels of the SpinGroup to those
 *  channels with the requested particle pair. This will allow a user to
 *  calculate observables of the compound system associated to different
 *  incident particle pairs.
 *
 *  @param[in] incident   the new incident particle pair
 */
void switchIncidentPair( const ParticlePair& incident ) {

  if ( incident.pairID() != this->incidentPair().pairID() ) {

    this->incident_ =
      determineIncidentChannels( incident.pairID(), this->channels_ );
    this->reactions_ =
      makeReactionIdentifiers( this->channels_, this->incident_.front() );

    //! @todo check for empty incident channels
  }
}

