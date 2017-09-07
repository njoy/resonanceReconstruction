class Lvalue : public lvalue::Type {
  using Parent = lvalue::Type;
  
  Lvalue() = delete;
  Lvalue( const Lvalue& ) = delete;
  Lvalue( Lvalue&& ) = delete;

  #include "resonanceReconstruction/breitWigner/multilevel/Lvalue/src/resonances.hpp"
  #include "resonanceReconstruction/breitWigner/multilevel/Lvalue/src/D.hpp"
  #include "resonanceReconstruction/breitWigner/multilevel/Lvalue/src/evaluate.hpp"

public:
  #include "resonanceReconstruction/breitWigner/multilevel/Lvalue/src/call.hpp"
};
