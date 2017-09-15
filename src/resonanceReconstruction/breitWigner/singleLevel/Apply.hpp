class Apply {
protected:
  template< template BreitWigner >
  static int nucleonNumber( const BreitWigner& slbw ){ 
    constexpr double neutronAMU = 1.00866491588;
    return std::lround( slbw.LStates().front().AWRI() * neutronAMU );
  }

  #include "resonanceReconstruction/breitWigner/singleLevel/Apply/src/competitiveWidth.hpp"
  #include "resonanceReconstruction/breitWigner/singleLevel/Apply/src/resonance.hpp"
  #include "resonanceReconstruction/breitWigner/singleLevel/Apply/src/lvalue.hpp"
  #include "resonanceReconstruction/breitWigner/singleLevel/Apply/src/build.hpp"

public:
  #include "resonanceReconstruction/breitWigner/singleLevel/Apply/src/call.hpp"  
};
