template< typename PsiChi, typename CompetitiveWidth = ZeroWidth >
auto evaluate( const PsiChi& kernel,
               const double channelRatio,
               const double scatteringRatio,
               const CompetitiveWidth& competitiveWidth = {} ) const {
  const auto penetrationShift = this->penetrationShift( channelRatio );
  const auto& penetrationFactor = penetrationShift[0];
  const auto& shiftFactor = penetrationShift[1];

  const auto phaseShift = this->phaseShift( scatteringRatio );

  const auto sinSquaredPhi = [&]{
    const auto sin = std::sin( phaseShift );
    return sin * sin;
  }();

  const auto sin2Phi = std::sin( 2. * phaseShift );

  const auto crossSections =
    this->resonances()
    | ranges::views::transform( [&]( const auto& resonance ){              
        return resonance( penetrationFactor,
                          shiftFactor,
                          sinSquaredPhi,
                          sin2Phi,
                          competitiveWidth( resonance ),
                          kernel ); } );

  const auto potentialScattering =
    pack( sinSquaredPhi * ( 2. * orbitalAngularMomentum + 1. ), 0.0, 0.0 );

  return ranges::accumulate( crossSections, potentialScattering );
}
