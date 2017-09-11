template< typename ChannelRadius, typename WaveNumber >
auto evaluate( const Quantity< ElectronVolts > energy,
               const ResonanceShape& kernel,
               const Quantity< RootBarn > channelRadius,
               const Quantity< RootBarn > scatteringRadius,
               ChannelRadius&& a ) const {
  const auto k = this->neutronWaveNumber();
  const Quantity< InvRootBarn > waveNumber = k( energy );
  const double channelRatio = channelRadius * waveNumber;
  const double scatteringRatio = scatteringRadius * waveNumber;

  const auto rho = [&]( Quantity< ElectronVolts > energy ){
    return k( energy ) * a( energy );
  };
  
  const auto crossSections =
    this->derived().lvalues()
    | ranges::view::transform
      ( [&]( auto&& lvalue ){ return lvalue( energy,
                                             kernel,
                                             channelRatio,
                                             scatteringRatio,
                                             channelRadius,
                                             this->targetSpin,
                                             rho ); } );
  
  const auto scaling = 4.0 * pi / ( waveNumber * waveNumber );
  const auto reduction = ranges::accumulate( crossSections,
                                             pack( 0.0, 0.0, 0.0 ) ).data; 

  return CrossSection( std::get<0>( reduction ) * scaling,
                       std::get<1>( reduction ) * scaling,
                       std::get<2>( reduction ) * scaling );
}
