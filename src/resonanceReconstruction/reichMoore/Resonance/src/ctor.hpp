Resonance( Quantity< ElectronVolts > energy,
           Quantity< ElectronVolts > neutronWidth,
           Quantity< ElectronVolts > captureWidth,
           Quantity< ElectronVolts > fissionWidthA,
           Quantity< ElectronVolts > fissionWidthB,
           double penetrationFactor,
           double flaggedStatisticalFactor ) :
  energy( energy ),
  rootNeutronWidth( sqrt( neutronWidth ) ),
  captureWidth( captureWidth ),
  rootFissionWidthA( (fissionWidthA.value < 0.0 ? -1 : 1)
                     * sqrt( std::abs( fissionWidthA ) ) ),
  rootFissionWidthB( (fissionWidthB.value < 0.0 ? -1 : 1)
                     * sqrt( std::abs( fissionWidthB ) ) ),
  inverseRootPenetrationFactor( 1. / std::sqrt( penetrationFactor ) ),
  flaggedStatisticalFactor( flaggedStatisticalFactor ){}
