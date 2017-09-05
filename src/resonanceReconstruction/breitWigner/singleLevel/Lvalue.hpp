class Lvalue : public lvalue::Type {
  using Parent = lvalue::Type;
  
  Lvalue() = delete;
  Lvalue( const Lvalue& ) = delete;
  Lvalue( Lvalue&& ) = delete;

  #include "resonanceReconstruction/breitWigner/singleLevel/Lvalue/src/resonances.hpp"
  #include "resonanceReconstruction/breitWigner/singleLevel/Lvalue/src/competitiveWidth.hpp"
  #include "resonanceReconstruction/breitWigner/singleLevel/Lvalue/src/evaluate.hpp"

public:
  #include "resonanceReconstruction/breitWigner/singleLevel/Lvalue/src/call.hpp"
};
