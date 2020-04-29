class Resonance : public resonance::Type {
  Resonance() = delete;
  Resonance( const Resonance& ) = delete;
  Resonance( Resonance&& ) = delete;
public:
  #include "resonanceReconstruction/breitWigner/multiLevel/Resonance/src/call.hpp"
};
