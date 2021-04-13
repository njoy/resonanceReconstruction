template< typename CompetitiveWidth = ZeroWidth >
auto evaluate( const ResonanceShape& kernel,
               const double channelRatio,
               const double scatteringRatio,
               const double targetSpin,
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

  const auto Jgroups =
    this->resonances()
    | ranges::views::group_by( []( const auto& reference, const auto& trial ){
        return reference.statisticalFactor == trial.statisticalFactor;
      } );

  const auto zero = pack( pack( 0., 0. ), 0., 0. );

  const auto crossSections =
    Jgroups
    | ranges::views::transform( [&]( const auto& Jgroup ){
        const auto Jresonances =
          Jgroup
          | ranges::views::transform(
              [&]( const auto& resonance ){
                return resonance( penetrationFactor,
                                  shiftFactor,
                                  competitiveWidth( resonance ),
                                  kernel ); } );

        const auto Jsum = ranges::accumulate( Jresonances, zero ).data;

        const auto& scatteringComponents = std::get<0>( Jsum ).data;
        const auto& first = std::get<0>( scatteringComponents );
        const auto& second = std::get<1>( scatteringComponents );

        const auto scattering = first * ( first - 2 * sinSquaredPhi )
                                + second * ( second + sin2Phi );
        const auto capture = std::get<1>( Jsum );
        const auto fission = std::get<2>( Jsum );

        return pack( scattering, capture, fission )
               * Jgroup.front().statisticalFactor;
      } );

  const double scatteringBase =
    this->D( targetSpin ) * sinSquaredPhi
    + this->statisticalFactorSum
      * ( sinSquaredPhi * sinSquaredPhi + 0.25 * sin2Phi * sin2Phi );

  const auto potentialScattering = pack( scatteringBase, 0., 0. );

  return ranges::accumulate( crossSections, potentialScattering );
}
