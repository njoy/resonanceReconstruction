struct Lvalue {
  
  #include "resonanceReconstruction/singleLevelBreitWigner/Base/Lvalue/Resonance.hpp"
  std::vector< Resonance > resonances;
  
  std::optional< Quantity<ElectronVolts > > weightedQValue;
  double target2CompoundWeightRatio;
  int orbitalAngularMomentum;

  #include "resonanceReconstruction/singleLevelBreitWigner/Base/Lvalue/src/waveNumber.hpp"
  #include "resonanceReconstruction/singleLevelBreitWigner/Base/Lvalue/src/penetrationShift.hpp"
  #include "resonanceReconstruction/singleLevelBreitWigner/Base/Lvalue/src/phaseShift.hpp"
  #include "resonanceReconstruction/singleLevelBreitWigner/Base/Lvalue/src/withCompetitiveWidth.hpp"
  #include "resonanceReconstruction/singleLevelBreitWigner/Base/Lvalue/src/withoutCompetitiveWidth.hpp"
  #include "resonanceReconstruction/singleLevelBreitWigner/Base/Lvalue/src/call.hpp"
};
