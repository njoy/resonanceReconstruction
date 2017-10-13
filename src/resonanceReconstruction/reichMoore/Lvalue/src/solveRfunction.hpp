static auto solveRfunction( const Matrix3x3& rMatrix,
                            const double cos2Phi,
                            const double sin2Phi,
                            const double statisticalFactor,
                            const double phaseShift ){
  const auto rMatrixEntry = 1. + rMatrix( 0, 0 );

  const double conjugateProduct =
    std::real( rMatrixEntry * std::conj( rMatrixEntry ) );

  const double denominator = 1. / conjugateProduct;

  const double real = std::real( rMatrixEntry ) * denominator;
  const double imaginary = -std::imag( rMatrixEntry ) * denominator;

  const std::complex<double> collisionMatrixElement =
    { cos2Phi * ( 2 * real - 1. ) + 2 * imaginary * sin2Phi,
      -sin2Phi * ( 2 * real - 1. ) + 2 * imaginary * cos2Phi };

  if ( ( std::abs( std::real( rMatrixEntry ) ) > 3e-4 )
       or ( std::abs( phaseShift ) > 3e-4 ) ){

    const auto conjugateComplement = 1. - std::conj( collisionMatrixElement );

    const double total =
      2 * statisticalFactor * std::real( conjugateComplement );

    const double elastic =
      statisticalFactor
      * ( std::real( conjugateComplement * std::conj( conjugateComplement ) ) );

    const double fission = 0.0;

    return pack( total, elastic, fission, statisticalFactor );
  } else {
    const auto expansion =
      2 * denominator * ( std::real( rMatrixEntry )
                          + conjugateProduct
                          + sin2Phi * std::imag( rMatrixEntry )
                          + phaseShift * phaseShift
                          * ( 1. - conjugateProduct ) );

    const double total = 2 * statisticalFactor * expansion;

    const double elastic =
      statisticalFactor * ( expansion * expansion
                            + std::imag( collisionMatrixElement )
                            * std::imag( collisionMatrixElement ) );

    const double fission = 0.0;

    return pack( total, elastic, fission, statisticalFactor );
  }
}
