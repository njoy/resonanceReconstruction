template< typename PsiChi, typename CompetitiveWidth >
CrossSection operator()( const Quantity<ElectronVolts> energy,
                         const PsiChi& kernel,
                         const Quantity<InvRootBarn> waveNumber,
                         const double channelRatio,
                         const double scatteringRatio,
                         const CompetitiveWidth& competitiveWidth ) const {
  const auto ps = this->penetrationShift( channelRatio );
  const auto& penetrationFactor = ps[0];
  const auto& shiftFactor = ps[1];

  const auto phaseShift = this->phaseShift( scatteringRatio );

  const auto trig = [&]() -> std::array<double, 3> {
    const auto sin = std::sin( phaseShift );
    const auto sinSquared = sin * sin;
    const auto cos = std::sqrt( 1. - sinSquared );
    const auto sin2 = 2. * sin * cos;
    const auto cos2 = 1. - 2. * sinSquared;
    return {{ sinSquared, sin2, cos2 }};
  }();
  
  const auto crossSections = 
    this->resonances
    | ranges::view::transform( [&]( auto&& resonance ){              
        return resonance( energy,
                          waveNumber,
                          penetrationFactor,
                          shiftFactor,
                          trig[1],
                          trig[2],
                          competitiveWidth( resonance ),
                          kernel ); } );

  const CrossSection potentialScattering =
    { 4. * pi * trig[0]
      * ( 2 * orbitalAngularMomentum + 1 )
      / ( waveNumber * waveNumber ),
      0.0 * barns,
      0.0 * barns };
  
  return ranges::accumulate( crossSections, potentialScattering );
} 

template< typename PsiChi,
          typename ChannelRadius,
          typename ScatteringRadius >
CrossSection operator()( const Quantity<ElectronVolts> energy,
                         const PsiChi& kernel,
                         ChannelRadius&& channelRadius,
                         ScatteringRadius&& scatteringRadius ) const {
  const auto waveNumber = this->waveNumber( energy );
  const auto channelRatio = waveNumber * channelRadius( energy );
  const auto scatteringRatio = waveNumber * scatteringRadius( energy );
  return ( this->weightedQValue ) ?
    (*this)( energy, kernel, waveNumber, channelRatio, scatteringRatio,
             this->withCompetitiveWidth( energy, channelRadius ) ) :
    (*this)( energy, kernel, waveNumber, channelRatio, scatteringRatio,
             this->withoutCompetitiveWidth() );
}

template< typename PsiChi, typename Radius >
CrossSection operator()( const Quantity<ElectronVolts> energy,
                         const PsiChi& kernel,
                         Radius&& radius ) const {
  const auto waveNumber = this->waveNumber( energy );
  const auto radius_ = radius( energy );
  const auto channelRatio = waveNumber * radius_;
  const auto scatteringRatio = waveNumber * radius_;
  return ( this->weightedQValue ) ?
    (*this)( energy, kernel, waveNumber, channelRatio, scatteringRatio,
             this->withCompetitiveWidth( energy, radius ) ) :
    (*this)( energy, kernel, waveNumber, channelRatio, scatteringRatio,
             this->withoutCompetitiveWidth() );
}
