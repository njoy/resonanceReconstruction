auto operator()( const Quantity<ElectronVolts> energy,
                 const double channelRatio,
                 const double scatteringRatio ) const {
  const double rootPenetrationFactor =
    std::sqrt( this->penetrationShift( channelRatio )[0] );

  const auto phaseShift = this->phaseShift( scatteringRatio );

  const auto sinSquaredPhi = [&]{
    const auto sin = std::sin( phaseShift );
    return sin * sin;
  }();

  const auto sin2Phi = std::sin( 2. * phaseShift );
  const auto cos2Phi = std::cos( 2. * phaseShift );

  const auto crossSections =
    this->resonances
    | ranges::view::group_by( []( auto&& left, auto&& right ){
        return left.flaggedStatisticalFactor == right.flaggedStatisticalFactor; } )
    | ranges::view::transform( [&]( const auto& group ){
        Matrix3x3 rMatrix = Matrix3x3::Zero();

        const bool hasFission =
          ranges::accumulate
          ( group | ranges::view::transform( [&]( const auto& resonance ) {
              return resonance( energy, rootPenetrationFactor, rMatrix );
            } ), false, std::logical_or<>{} );

        const double statisticalFactor =
          std::abs( group.front().flaggedStatisticalFactor );

        return ( hasFission ) ?
          solveRmatrix( rMatrix, cos2Phi, sin2Phi, statisticalFactor ) :
          solveRfunction( rMatrix, cos2Phi, sin2Phi, statisticalFactor, phaseShift );
      } );

  const auto sum =
    ranges::accumulate( crossSections, pack( 0., 0., 0., 0. ) ).data;

  const auto unaccountedPotentialScattering =
    ( ( 2 * this->orbitalAngularMomentum + 1. ) - std::get<3>( sum ) )
    * 4 * sinSquaredPhi;

  const auto elastic = std::get<1>(sum) + unaccountedPotentialScattering;
  const auto capture = std::get<0>(sum) - std::get<1>(sum) - std::get<2>(sum);
  const auto fission = std::get<2>(sum);
  return pack( elastic, capture, fission );
}
