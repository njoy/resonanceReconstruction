bool operator()( const Quantity< ElectronVolts > energy,
                 const double rootPenetrationFactor,
                 Matrix3x3& K ) const {
  const auto rootNeutronWidth = this->rootNeutronWidth
                                * this->inverseRootPenetrationFactor
                                * rootPenetrationFactor;

  const auto energyDifference = (this->energy - energy);
  const auto squareEnergyDifference = energyDifference * energyDifference;
  const auto squareCaptureWidth = this->captureWidth * this->captureWidth;

  const auto denominator =
    1. / ( 4 * squareEnergyDifference + squareCaptureWidth );

  const auto realFactor = this->captureWidth * denominator;
  const auto imaginaryFactor = -2 * energyDifference * denominator;

  auto accumulate = [&]( auto& element, const auto width0, const auto width1 ){
    const auto product = width0 * width1;
    const double real = product * realFactor;
    const double imaginary = product * imaginaryFactor;
    element += std::complex<double>{real, imaginary};
  };

  accumulate( K(0,0), rootNeutronWidth, rootNeutronWidth );

  if ( this->rootFissionWidthA.value != 0.0
       or this->rootFissionWidthB.value != 0.0 ){
    accumulate( K( 0, 1 ), rootNeutronWidth, this->rootFissionWidthA );
    accumulate( K( 0, 2 ), rootNeutronWidth, this->rootFissionWidthB );
    accumulate( K( 1, 1 ), this->rootFissionWidthA, this->rootFissionWidthA );
    accumulate( K( 1, 2 ), this->rootFissionWidthA, this->rootFissionWidthB );
    accumulate( K( 2, 2 ), this->rootFissionWidthB, this->rootFissionWidthB );
    return true;
  }
  return false;
}
