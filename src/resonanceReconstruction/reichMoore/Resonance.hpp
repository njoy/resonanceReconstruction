class Resonance {
public:
  Quantity< ElectronVolts > energy;
  Quantity< RootElectronVolts > rootNeutronWidth;
  Quantity< ElectronVolts > captureWidth;
  Quantity< RootElectronVolts > rootFissionWidthA;
  Quantity< RootElectronVolts > rootFissionWidthB;
  double inverseRootPenetrationFactor;
  double flaggedStatisticalFactor;

protected:
  bool padding = true; // increase struct size to 64 bytes for optimal alignment

public:
  #include "resonanceReconstruction/reichMoore/Resonance/src/ctor.hpp"
  #include "resonanceReconstruction/reichMoore/Resonance/src/call.hpp"

};
