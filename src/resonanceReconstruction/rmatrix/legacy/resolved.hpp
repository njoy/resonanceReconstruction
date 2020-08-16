namespace resolved {

  // data collections for the legacy resolved resonance reconstruction
  using Components = Data< double >;

  // resolved resonance data
  #include "resonanceReconstruction/rmatrix/legacy/resolved/Resonance.hpp"
  using ResonanceTable = ResonanceTableBase< Resonance >;
}
