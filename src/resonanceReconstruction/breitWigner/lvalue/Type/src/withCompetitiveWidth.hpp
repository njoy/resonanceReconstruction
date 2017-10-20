template< typename ChannelRadius >
auto withCompetitiveWidth( Quantity<ElectronVolts> energy,
                           ChannelRadius&& channelRadius ) const {
  energy += *(this->weightedQValue);
  const auto waveNumber = this->waveNumber( energy );
  const auto channelRatio = waveNumber * channelRadius( energy );
  const auto penetrationFactor = this->penetrationShift( channelRatio )[0];

  return [ penetrationFactor ]( auto&& resonance ){
    return resonance.weightedCompetitiveWidth * penetrationFactor;
  };
}
