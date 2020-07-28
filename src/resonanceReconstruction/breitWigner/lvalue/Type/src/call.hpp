template< typename PsiChi, typename CompetitiveWidth >
auto operator()( const PsiChi& kernel,
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
    | ranges::view::transform( [&]( const auto& resonance ){
        return resonance( penetrationFactor,
                          shiftFactor,
                          trig[1],
                          trig[2],
                          competitiveWidth( resonance ),
                          kernel ); } );

  const auto potentialScattering =
    pack( trig[0] * ( 2 * orbitalAngularMomentum + 1 ), 0.0, 0.0 );

  return ranges::accumulate( crossSections, potentialScattering );
}

template< typename PsiChi, typename ChannelRadius >
auto operator()( const Quantity<ElectronVolts> energy,
                 const PsiChi& kernel,
                 const double channelRatio,
                 const double scatteringRatio,
                 ChannelRadius&& channelRadius ) const {
  return ( this->weightedQValue ) ?
    (*this)( kernel, channelRatio, scatteringRatio,
             this->withCompetitiveWidth( energy, channelRadius ) ) :
    (*this)( kernel, channelRatio, scatteringRatio,
             this->withoutCompetitiveWidth() );
}
