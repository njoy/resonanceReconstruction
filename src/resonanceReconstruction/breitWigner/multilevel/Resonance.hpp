class Resonance : public resonance::Type {
  Resonance() = delete;
  Resonance( const Resonance& ) = delete;
  Resonance( Resonance&& ) = delete;
public:
  #include "resonanceReconstruction/breitWigner/multilevel/Resonance/src/call.hpp"
};
