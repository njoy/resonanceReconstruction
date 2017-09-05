template< typename PsiChi, typename CompetitiveWidth = ZeroWidth >
double evaluate( const PsiChi& kernel,
                 const double channelRatio,
                 const double scatteringRatio,
                 const CompetitiveWidth& competitiveWidth = {} ) const {
  const auto penetrationShift = this->penetrationShift( channelRatio );
  const auto& penetrationFactor = penetrationShift[0];
  const auto& shiftFactor = penetrationShift[1];

  const auto phaseShift = this->phaseShift( scatteringRatio );

  const auto trig = [&]() -> std::array< double, 3 > {
    const auto sin = std::sin( phaseShift );
    const auto sinSquared = sin * sin;
    const auto sin2 = std::sin( 2. * phaseShift );
    const auto cos2 = std::cos( 2. * phaseShift );
    return {{ sinSquared, sin2, cos2 }};
  }();

  const auto crossSections = 
    this->resonances()
    | ranges::view::transform( [&]( const auto& resonance ){              
        return resonance( penetrationFactor,
                          shiftFactor,
                          trig[1],
                          trig[2],
                          competitiveWidth( resonance ),
                          kernel ); } );

  const double potentialScattering =
    pack( trig[0] * ( 2. * orbitalAngularMomentum + 1. ), 0.0, 0.0 );

  return ranges::accumulate( crossSections, potentialScattering );
} 
