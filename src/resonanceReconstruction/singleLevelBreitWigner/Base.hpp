struct Base {
  #include "resonanceReconstruction/singleLevelBreitWigner/Base/CrossSection.hpp"
  #include "resonanceReconstruction/singleLevelBreitWigner/Base/Lvalue.hpp"

  struct EnergyRange {
    Quantity< ElectronVolts > lowerLimit;
    Quantity< ElectronVolts > upperLimit;
  };

  std::vector< Lvalue > lvalues;
  EnergyRange energyRange;
  double nucleonNumber;

  Base( std::vector< Lvalue >&& lvalues,
        EnergyRange energyRange,
        int nucleonNumber ) :
    lvalues( std::move(lvalues) ),
    energyRange( energyRange ),
    nucleonNumber( nucleonNumber ){}

  #include "resonanceReconstruction/singleLevelBreitWigner/Base/src/psiChi.hpp"
  #include "resonanceReconstruction/singleLevelBreitWigner/Base/src/evaluate.hpp"
  #include "resonanceReconstruction/singleLevelBreitWigner/Base/src/isTemperature.hpp"
  #include "resonanceReconstruction/singleLevelBreitWigner/Base/src/call.hpp"
};
