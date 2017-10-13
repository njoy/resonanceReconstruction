class Lvalue {
public:
  using AP = decltype( radius(1.0) );

protected:
  std::optional< AP > ap_;
  std::vector< Resonance > resonances;
  int orbitalAngularMomentum;

public:
  #include "resonanceReconstruction/reichMoore/Lvalue/src/ctor.hpp"
  #include "resonanceReconstruction/reichMoore/Lvalue/src/ap.hpp"
  #include "resonanceReconstruction/reichMoore/Lvalue/src/penetrationShift.hpp"
  #include "resonanceReconstruction/reichMoore/Lvalue/src/phaseShift.hpp"
  #include "resonanceReconstruction/reichMoore/Lvalue/src/solveRMatrix.hpp"
  #include "resonanceReconstruction/reichMoore/Lvalue/src/solveRfunction.hpp"
  #include "resonanceReconstruction/reichMoore/Lvalue/src/call.hpp"
};
