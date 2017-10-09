class Type {
protected:
  std::vector< lvalue::Type > lvalues;
  double target2CompoundWeightRatio;
  EnergyRange energyRange;
  double nucleonNumber;
  
  #include "resonanceReconstruction/breitWigner/Type/src/ctor.hpp"  
  #include "resonanceReconstruction/breitWigner/Type/src/neutronWaveNumber.hpp"
};
