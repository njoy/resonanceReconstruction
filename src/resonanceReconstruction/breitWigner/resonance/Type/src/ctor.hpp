Type( const Quantity< ElectronVolts > energy,
      const Quantity< ElectronVolts > neutronWidth,
      const Quantity< ElectronVolts > captureWidth,
      const Quantity< ElectronVolts > fissionWidth,
      const Quantity< ElectronVolts > weightedCompetitiveWidth,
      const double penetrationFactor,
      const double shiftFactor,
      const double statisticalFactor ) :
  energy( energy ),
  neutronWidth( neutronWidth ),
  captureWidth( captureWidth ),
  fissionWidth( fissionWidth ),
  weightedCompetitiveWidth( weightedCompetitiveWidth ),
  inversePenetrationFactor( 1. / penetrationFactor ),
  shiftFactor( shiftFactor ),
  statisticalFactor( statisticalFactor ){}
