namespace unresolved {

  // data collections for the legacy unresolved resonance reconstruction
  using Widths = Data< Width >;
  using FluctuationIntegrals = Data< FluctuationIntegral >;
  using Degrees = Data< unsigned int >;

  // unresolved resonance data
  #include "resonanceReconstruction/rmatrix/legacy/unresolved/Resonance.hpp"
  #include "resonanceReconstruction/rmatrix/legacy/unresolved/ResonanceTable.hpp"

  // function to calculate the fluctuation integrals
  #include "resonanceReconstruction/rmatrix/legacy/unresolved/src/calculateFluctuationIntegrals.hpp"

  // unresolved spin group and compound system
  #include "resonanceReconstruction/rmatrix/legacy/unresolved/SpinGroup.hpp"
  using CompoundSystem = CompoundSystemBase< unresolved::SpinGroup >;
}
