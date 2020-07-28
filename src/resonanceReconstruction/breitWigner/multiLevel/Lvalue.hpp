class Lvalue : public lvalue::Type {
  using Parent = lvalue::Type;

  Lvalue() = delete;
  Lvalue( const Lvalue& ) = delete;
  Lvalue( Lvalue&& ) = delete;

  #include "resonanceReconstruction/breitWigner/multiLevel/Lvalue/src/resonances.hpp"
  #include "resonanceReconstruction/breitWigner/multiLevel/Lvalue/src/D.hpp"
  #include "resonanceReconstruction/breitWigner/multiLevel/Lvalue/src/evaluate.hpp"

public:
  #include "resonanceReconstruction/breitWigner/multiLevel/Lvalue/src/call.hpp"
};
