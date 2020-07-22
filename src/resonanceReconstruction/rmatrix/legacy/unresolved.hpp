namespace unresolved {

  // data collections for the legacy unresolved resonance reconstruction
  #include "resonanceReconstruction/rmatrix/legacy/unresolved/Data.hpp"
  using Widths = Data< Width >;
  using FluctuationIntegrals = Data< FluctuationIntegral >;
  using CrossSections = Data< CrossSection >;
  using Degrees = Data< unsigned int >;

  // function to calculate the fluctuation integrals
  #include "resonanceReconstruction/rmatrix/legacy/unresolved/src/calculateFluctuationIntegrals.hpp"
}
