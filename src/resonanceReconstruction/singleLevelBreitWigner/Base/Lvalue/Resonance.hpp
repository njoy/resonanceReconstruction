struct Resonance {
  Quantity<ElectronVolts> energy;
  Quantity<ElectronVolts> neutronWidth;
  Quantity<ElectronVolts> captureWidth;
  Quantity<ElectronVolts> fissionWidth;
  Quantity<ElectronVolts> weightedCompetitiveWidth;
  double inversePenetrationFactor;
  double shiftFactor;
  double statisticalFactor;

  #include "resonanceReconstruction/singleLevelBreitWigner/Base/Lvalue/Resonance/src/call.hpp"
};
