namespace resolved {

  // resolved resonance data
  #include "resonanceReconstruction/rmatrix/legacy/resolved/Resonance.hpp"
  using ResonanceTable = ResonanceTableBase< resolved::Resonance >;

  // legacy resolved formalism options
  struct SingleLevelBreitWigner {};
  struct MultiLevelBreitWigner {};

  // resolved spin group and compound system
  #include "resonanceReconstruction/rmatrix/legacy/resolved/SpinGroup.hpp"
  template < typename Formalism > using CompoundSystem = CompoundSystemBase< resolved::SpinGroup< Formalism > >;
}
