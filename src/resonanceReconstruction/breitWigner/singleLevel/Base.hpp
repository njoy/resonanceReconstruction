class Base : public breitWigner::Type {
protected:
  using Parent = breitWigner::Type;
  
  #include "resonanceReconstruction/breitWigner/singleLevel/Base/src/psiChi.hpp"
  #include "resonanceReconstruction/breitWigner/singleLevel/Base/src/lvalues.hpp"

  using breitWigner::Type::Type;
};
