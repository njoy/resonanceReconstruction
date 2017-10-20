static double
statisticalFactorSummation( const double targetSpin,
                            const double orbitalAngularMomentum ){
  const double jmin =
    std::abs( std::abs( targetSpin - orbitalAngularMomentum ) - 0.5 );
  const double jmax = targetSpin + orbitalAngularMomentum + 0.5;

  return ( jmax - jmin + 1. ) * ( jmax + jmin + 1. ) / ( 4 * targetSpin + 2. );
}
