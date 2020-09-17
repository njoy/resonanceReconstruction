namespace resolved {

  // resolved resonance data
  #include "resonanceReconstruction/rmatrix/legacy/resolved/Resonance.hpp"
  using ResonanceTable = ResonanceTableBase< resolved::Resonance >;

  // legacy resolved formalism options
  struct SingleLevelBreitWigner {};
  struct MultiLevelBreitWigner {};

  // resolved spin group and compound system
  #include "resonanceReconstruction/rmatrix/legacy/resolved/SpinGroup.hpp"
  #include "resonanceReconstruction/rmatrix/legacy/resolved/CompoundSystem.hpp"
}
