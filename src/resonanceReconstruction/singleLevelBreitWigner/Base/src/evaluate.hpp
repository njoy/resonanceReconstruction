template< typename PsiChi, typename WaveNumber, typename ChannelRadius >
auto evaluate( Quantity<ElectronVolts> energy,
               const PsiChi& kernel,
               const WaveNumber waveNumber,
               const double channelRatio,
               const double scatteringRatio,
               ChannelRadius&& channelRadius ) const {  
  const auto crossSections =
    this->lvalues
    | ranges::view::transform
      ( [&]( auto&& lvalue ){
        return lvalue( energy,
                       kernel,
                       channelRatio,
                       scatteringRatio,
                       channelRadius ); } );
  
  const auto scaling = 4 * pi / ( waveNumber * waveNumber ); 
  const auto reduction = ranges::accumulate( crossSections,
                                             pack( 0.0, 0.0, 0.0 ) ); 

  return CrossSection( std::get<0>( reduction.data ) * scaling,
                       std::get<1>( reduction.data ) * scaling,
                       std::get<2>( reduction.data ) * scaling );
}

template< typename PsiChi, typename Radius >
auto evaluate( Quantity<ElectronVolts> energy,
               const PsiChi& kernel,
               Radius&& radius ) const {
  /* we've assumed the atomic weight ratio is consistent over lvalues */
  const auto waveNumber = this->lvalues.front().waveNumber( energy );
  const auto channelRatio = waveNumber * radius( energy );
  const auto scatteringRatio = waveNumber * radius( energy );
  return this->evaluate( energy,
                         kernel,
                         waveNumber,
                         channelRatio,
                         scatteringRatio,
                         radius );
}

template< typename PsiChi, typename ChannelRadius, typename ScatteringRadius >
auto evaluate( Quantity<ElectronVolts> energy,
               // from child call
               const PsiChi& kernel,
               ChannelRadius&& channelRadius,
               ScatteringRadius&& scatteringRadius ) const {
  /* we've assumed the atomic weight ratio is consistent over lvalues */
  const auto waveNumber = this->lvalues.front().waveNumber( energy );
  const auto channelRatio = waveNumber * channelRadius( energy );
  const auto scatteringRatio = waveNumber * scatteringRadius( energy );
  return this->evaluate( energy,
                         kernel,
                         waveNumber,
                         channelRatio,
                         scatteringRatio,
                         channelRadius );
}
