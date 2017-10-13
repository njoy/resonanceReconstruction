static auto solveRmatrix( Matrix3x3& rMatrix,
                          const double cos2Phi,
                          const double sin2Phi,
                          const double statisticalFactor ){
  rMatrix.triangularView< Eigen::Lower >() = rMatrix.transpose();

  const Matrix3x3 collisionMatrix =
    2 * ( Matrix3x3::Identity() + rMatrix ).inverse() - Matrix3x3::Identity();

  const auto elasticTerm =
    std::complex<double>{ cos2Phi, sin2Phi } * collisionMatrix( 0, 0 );

  const double total =
    2 * statisticalFactor * ( 1. - std::real( elasticTerm ) );

  const double elastic =
    statisticalFactor * ( std::pow( 1. - std::real( elasticTerm ), 2 )
                          + std::pow( std::imag( elasticTerm ), 2 ) );

  const double fission =
    4 * statisticalFactor
    * ( std::real( collisionMatrix(0,1)
                   * std::conj( collisionMatrix(0,1) ) )
        + std::real( collisionMatrix(0,2)
                     * std::conj( collisionMatrix(0,2) ) ) );

  return pack( total, elastic, fission, statisticalFactor );
}
