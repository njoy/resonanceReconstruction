namespace rmatrix {

  // quantum numbers
  using OrbitalAngularMomentum = unsigned int;
  using Spin = double;
  using TotalAngularMomentum = double;
  using Parity = short;

  // auxiliary functions for the quantum numbers
  #include "resonanceReconstruction/rmatrix/src/possibleChannelSpinValues.hpp"
  #include "resonanceReconstruction/rmatrix/src/possibleChannelTotalAngularMomentumValues.hpp"

  // particle and channel types, functions dependent on those types
  struct Neutron {};
  struct Photon {};
  struct ChargedParticle {};
  struct Fission {};
  #include "resonanceReconstruction/rmatrix/src/coh3-coulomb.hpp"
  #include "resonanceReconstruction/rmatrix/src/calculatePenetrability.hpp"
  #include "resonanceReconstruction/rmatrix/src/calculateShiftFactor.hpp"
  #include "resonanceReconstruction/rmatrix/src/calculatePhaseShift.hpp"
  #include "resonanceReconstruction/rmatrix/src/calculateCoulombPhaseShift.hpp"

  // identifiers
  using ParticleID = std::string;
  using ParticlePairID = std::string;
  using ReactionID = std::string;
  using ChannelID = std::string;

  // matrix
  template < typename T > using Matrix = Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic >;

  // R-Matrix boundary condition and options
  using BoundaryCondition = double;
  struct ShiftFactor {};
  struct Constant {};
  #include "resonanceReconstruction/rmatrix/src/calculateLValue.hpp"

  // R-matrix components (independent of formalism)
  #include "resonanceReconstruction/rmatrix/Particle.hpp"
  #include "resonanceReconstruction/rmatrix/ParticlePair.hpp"
  #include "resonanceReconstruction/rmatrix/ChannelQuantumNumbers.hpp"
  #include "resonanceReconstruction/rmatrix/ChannelRadii.hpp"
  #include "resonanceReconstruction/rmatrix/Channel.hpp"
  #include "resonanceReconstruction/rmatrix/Resonance.hpp"
  #include "resonanceReconstruction/rmatrix/ResonanceTable.hpp"

  // R-Matrix formalism options
  struct ReichMoore {};
  struct GeneralRMatrix {};

  // spin group and compound system
  #include "resonanceReconstruction/rmatrix/SpinGroup.hpp"
  #include "resonanceReconstruction/rmatrix/CompoundSystem.hpp"

  // make components from ENDF
  #include "resonanceReconstruction/rmatrix/src/makeParticlePairs.hpp"
  #include "resonanceReconstruction/rmatrix/src/makeParticleChannels.hpp"
  #include "resonanceReconstruction/rmatrix/src/makeResonanceTable.hpp"
  #include "resonanceReconstruction/rmatrix/src/makeSpinGroup.hpp"
  #include "resonanceReconstruction/rmatrix/src/makeCompoundSystem.hpp"

  // make the compound system variant from ENDF
  using CompoundSystemVariant =
      std::variant< CompoundSystem< ReichMoore, ShiftFactor >,
                    CompoundSystem< ReichMoore, Constant >/*,
                    CompoundSystem< GeneralRMatrix, ShiftFactor >,
                    CompoundSystem< GeneralRMatrix, Constant >*/ >;
  #include "resonanceReconstruction/rmatrix/Reconstructor.hpp"
  #include "resonanceReconstruction/rmatrix/src/fromENDF.hpp"
}
