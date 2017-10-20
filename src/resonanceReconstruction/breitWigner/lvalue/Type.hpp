class Type {
protected:
  std::vector< resonance::Type > resonances;
  std::optional< Quantity< ElectronVolts > > weightedQValue;
  int orbitalAngularMomentum;
  double statisticalFactorSum;

  #include "resonanceReconstruction/breitWigner/lvalue/Type/src/statisticalFactorSummation.hpp"

public:
  #include "resonanceReconstruction/breitWigner/lvalue/Type/src/ctor.hpp"
  #include "resonanceReconstruction/breitWigner/lvalue/Type/src/penetrationShift.hpp"
  #include "resonanceReconstruction/breitWigner/lvalue/Type/src/phaseShift.hpp"
  #include "resonanceReconstruction/breitWigner/lvalue/Type/src/competitiveWidth.hpp"
};
