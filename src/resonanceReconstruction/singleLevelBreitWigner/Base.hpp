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

  auto makePsiChi( Quantity<ElectronVolts> e ){
    return [e]( auto ep, auto gammaInv )->std::array<double, 2>{
      auto x = 2 * ( e - ep ) * gammaInv;
      auto psi = 1 / ( 1 + ( x * x ) );
      return std::array<double,2>{ psi, x * psi };
    };
  }

  Base( std::vector< Lvalue >&& lvalues,
        EnergyRange energyRange,
        int atomicNumber ) :
    lvalues( std::move(lvalues) ),
    energyRange( energyRange ),
    atomicNumber( atomicNumber ){}
  
  // #include "resonanceReconstruction/singleLevelBreitWigner/Base/src/call.hpp"
};
