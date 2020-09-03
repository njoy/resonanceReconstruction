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
  #include "resonanceReconstruction/rmatrix/src/horner.hpp"
  #include "resonanceReconstruction/rmatrix/src/coh3-coulomb.hpp"
  #include "resonanceReconstruction/rmatrix/src/calculatePenetrability.hpp"
  #include "resonanceReconstruction/rmatrix/src/calculateShiftFactor.hpp"
  #include "resonanceReconstruction/rmatrix/src/calculatePhaseShift.hpp"
  #include "resonanceReconstruction/rmatrix/src/calculateCoulombPhaseShift.hpp"

  // identifiers
  using ParticleID = elementary::ParticleID;
  using ParticlePairID = elementary::ParticlePairID;
  using ReactionID = elementary::ReactionID;
  using ChannelID = std::string;
  using ReactionChannelID = std::string;

  // matrix
  template < typename T > using Matrix = Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic >;
  template < typename T > using DiagonalMatrix = Eigen::DiagonalMatrix< T, Eigen::Dynamic >;

  // utility code
  #include "resonanceReconstruction/rmatrix/Table.hpp"
  #include "resonanceReconstruction/rmatrix/overload.hpp"

  // R-Matrix boundary condition and options
  using BoundaryCondition = double;
  struct ShiftFactor {};
  struct Constant {};
  #include "resonanceReconstruction/rmatrix/LMatrixCalculator.hpp"

  // R-matrix components (independent of formalism)
  #include "resonanceReconstruction/rmatrix/Particle.hpp"
  #include "resonanceReconstruction/rmatrix/ParticlePair.hpp"
  #include "resonanceReconstruction/rmatrix/ChannelQuantumNumbers.hpp"
  #include "resonanceReconstruction/rmatrix/ChannelRadiusTable.hpp"
  #include "resonanceReconstruction/rmatrix/ChannelRadii.hpp"
  #include "resonanceReconstruction/rmatrix/Channel.hpp"
  #include "resonanceReconstruction/rmatrix/ParticleChannel.hpp"
  #include "resonanceReconstruction/rmatrix/ParticleChannelData.hpp"

  // resolved resonance information
  #include "resonanceReconstruction/rmatrix/Resonance.hpp"
  #include "resonanceReconstruction/rmatrix/ResonanceTable.hpp"

  // legacy resonance reconstruction
  #include "resonanceReconstruction/rmatrix/legacy.hpp"

  // R-Matrix formalism options
  struct ReichMoore {};
  struct GeneralRMatrix {};
  #include "resonanceReconstruction/rmatrix/RLMatrixCalculator.hpp"

  // spin group and compound system
  #include "resonanceReconstruction/rmatrix/SpinGroup.hpp"
  #include "resonanceReconstruction/rmatrix/CompoundSystem.hpp"

  // make components from ENDF
  #include "resonanceReconstruction/rmatrix/src/makeQuantumNumbers.hpp"
  #include "resonanceReconstruction/rmatrix/src/makeParticlePairs.hpp"
  #include "resonanceReconstruction/rmatrix/src/makeParticleChannels.hpp"
  #include "resonanceReconstruction/rmatrix/src/makeParticleChannelData.hpp"
  #include "resonanceReconstruction/rmatrix/src/makeChannelRadiusTable.hpp"
  #include "resonanceReconstruction/rmatrix/src/makeChannelRadii.hpp"
  #include "resonanceReconstruction/rmatrix/src/makeCompoundSystem.hpp"
  #include "resonanceReconstruction/rmatrix/src/makeReichMooreChannelData.hpp"
  #include "resonanceReconstruction/rmatrix/src/makeReichMooreCompoundSystem.hpp"
  #include "resonanceReconstruction/rmatrix/src/makeLegacyUnresolvedResonanceTable.hpp"
  #include "resonanceReconstruction/rmatrix/src/makeLegacyUnresolvedSpinGroups.hpp"
  #include "resonanceReconstruction/rmatrix/src/makeLegacyUnresolvedCompoundSystem.hpp"

  // make the compound system variant from ENDF
  using CompoundSystemVariant =
      std::variant< CompoundSystem< ReichMoore, ShiftFactor >,
                    CompoundSystem< ReichMoore, Constant >,
                    //CompoundSystem< GeneralRMatrix, ShiftFactor >,
                    //CompoundSystem< GeneralRMatrix, Constant >,
                    legacy::unresolved::CompoundSystem >;
  #include "resonanceReconstruction/rmatrix/Reconstructor.hpp"
  #include "resonanceReconstruction/rmatrix/src/fromENDF.hpp"
}
