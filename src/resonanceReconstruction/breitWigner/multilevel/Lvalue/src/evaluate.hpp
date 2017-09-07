template< typename PsiChi, typename CompetitiveWidth = ZeroWidth >
auto evaluate( const PsiChi& kernel,
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

  const auto Jgroups =
    this->resonances()
    | ranges::view::group_by( []( const auto& reference, const auto& trial ){
        return reference.statisticalFactor == trial.statisticalFactor; } );

  const auto base = pack( pack( trig[2] - 1., trig[1] ), 0.0, 0.0 );
  
  const auto crossSections = 
    Jgroups
    | ranges::view::transform( [&]( const auto& Jgroup ){
        const auto Jresonances =
          Jgroup
          | ranges::view::transform(
              [&]( const auto& resonance ){
                return resonance( penetrationFactor,
                                  shiftFactor,
                                  competitiveWidth( resonance ),
                                  kernel ); } );
        
        const auto Jsum = ranges::accumulate( Jresonances, base ).data;
        const auto& scatteringComponents = std::get<0>( Jsum );
        const auto& first = std::get<0>( scatteringComponents );
        const auto& second = std::get<1>( scatteringComponents );
        
        const auto scattering = ( first * first + second * second );
        const auto capture = std::get<1>( Jsum );
        const auto fission = std::get<2>( Jsum );
        
        const auto statisticalFactor = Jvalue.front().statisticalFactor;
        return pack( scattering, capture, fission ) * statisticalFactor;
      } );

  const auto potentialScattering =
    pack( 2. * this->D() * (trig[2] - 1.), 0., 0. );
  
  return ranges::accumulate( crossSections, potentialScattering );
} 
