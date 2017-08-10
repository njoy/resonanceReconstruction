struct Base {
  #include "resonanceReconstruction/singleLevelBreitWigner/Base/CrossSection.hpp"
  #include "resonanceReconstruction/singleLevelBreitWigner/Base/Lvalue.hpp"

  struct EnergyRange {
    Quantity< ElectronVolts > lowerLimit;
    Quantity< ElectronVolts > upperLimit;
  };

  std::vector< Lvalue > lvalues;
  EnergyRange energyRange;
  double atomicNumber;

  Base( std::vector< Lvalue >&& lvalues,
        EnergyRange energyRange,
        int atomicNumber ) :
    lvalues( std::move(lvalues) ),
    energyRange( energyRange ),
    atomicNumber( atomicNumber ){}
  
  // #include "resonanceReconstruction/singleLevelBreitWigner/Base/src/call.hpp"
};
