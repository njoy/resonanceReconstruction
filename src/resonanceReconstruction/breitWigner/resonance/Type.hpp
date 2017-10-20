class Type {
public:
  Quantity< ElectronVolts > energy;
  Quantity< ElectronVolts > neutronWidth;
  Quantity< ElectronVolts > captureWidth;
  Quantity< ElectronVolts > fissionWidth;
  Quantity< ElectronVolts > weightedCompetitiveWidth;
  double inversePenetrationFactor;
  double shiftFactor;
  double statisticalFactor;

  #include "resonanceReconstruction/breitWigner/resonance/Type/src/ctor.hpp"
};
