static auto solveRmatrix( Matrix3x3& rMatrix,
                          const double cos2Phi,
                          const double sin2Phi,
                          const double statisticalFactor ){
  rMatrix.triangularView< Eigen::Lower >() = rMatrix.transpose();

  const Matrix3x3 collisionMatrix =
    std::complex<double>( cos2Phi, -sin2Phi )
    * ( 2. * ( Matrix3x3::Identity() + rMatrix ).inverse()
        - Matrix3x3::Identity() );

  const double total =
    2 * statisticalFactor * ( 1. - std::real( collisionMatrix(0,0) ) );

  const double elastic =
    statisticalFactor * std::pow( std::abs( 1. - collisionMatrix(0,0) ), 2 );

  const double fission =
    statisticalFactor
    * ( std::abs( std::pow( collisionMatrix(0,1), 2 ) )
        + std::abs( std::pow( collisionMatrix(0,2), 2 ) ) );

  return pack( total, elastic, fission, statisticalFactor );
}
