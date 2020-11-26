#include "resonanceReconstruction/quantities.hpp"
#include "resonanceReconstruction/rmatrix/options.hpp"

// utility code
#include "resonanceReconstruction/rmatrix/Table.hpp"
#include "utility/overload.hpp"

// R-matrix components (independent of formalism)
#include "resonanceReconstruction/rmatrix/QuantumNumbers.hpp"
#include "resonanceReconstruction/rmatrix/ChannelQuantumNumbers.hpp"
#include "resonanceReconstruction/rmatrix/Particle.hpp"
#include "resonanceReconstruction/rmatrix/ParticlePair.hpp"
#include "resonanceReconstruction/rmatrix/ChannelRadiusTable.hpp"
#include "resonanceReconstruction/rmatrix/ChannelRadii.hpp"

// auxiliary functions for the quantum numbers
#include "resonanceReconstruction/rmatrix/possibleChannelSpinValues.hpp"
#include "resonanceReconstruction/rmatrix/possibleChannelTotalAngularMomentumValues.hpp"

#ifndef NJOY_R2_RMATRIX
#define NJOY_R2_RMATRIX

// system includes

// other includes

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {

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
  using ReactionType = elementary::ReactionType;
  using ReactionID = elementary::ReactionID;
  using ChannelID = std::string;
  using ReactionChannelID = std::string;

  // matrix
  template < typename T > using Matrix = Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic >;
  template < typename T > using DiagonalMatrix = Eigen::DiagonalMatrix< T, Eigen::Dynamic >;

  // map
  template < typename Key, typename Value > using Map = std::map< Key, Value >;

  // R-Matrix boundary condition
  using BoundaryCondition = double;
  #include "resonanceReconstruction/rmatrix/LMatrixCalculator.hpp"

  // R-matrix components (independent of formalism)
  #include "resonanceReconstruction/rmatrix/Channel.hpp"
  #include "resonanceReconstruction/rmatrix/ParticleChannel.hpp"
  #include "resonanceReconstruction/rmatrix/ParticleChannelData.hpp"

  // resonance information
  #include "resonanceReconstruction/rmatrix/Resonance.hpp"
  #include "resonanceReconstruction/rmatrix/ResonanceTable.hpp"

  // formalism calculators
  #include "resonanceReconstruction/rmatrix/RLMatrixCalculator.hpp"

  // spin group and compound system
  #include "resonanceReconstruction/rmatrix/SpinGroup.hpp"
  #include "resonanceReconstruction/rmatrix/CompoundSystem.hpp"

  // legacy resonance reconstruction
  #include "resonanceReconstruction/rmatrix/legacy.hpp"

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
  #include "resonanceReconstruction/rmatrix/src/makeLegacyBreitWignerSpinGroups.hpp"
  #include "resonanceReconstruction/rmatrix/src/makeLegacyBreitWignerCompoundSystem.hpp"
  #include "resonanceReconstruction/rmatrix/src/makeLegacyUnresolvedResonanceTable.hpp"
  #include "resonanceReconstruction/rmatrix/src/makeLegacyUnresolvedSpinGroups.hpp"
  #include "resonanceReconstruction/rmatrix/src/makeLegacyUnresolvedCompoundSystem.hpp"

  // make the compound system variant from ENDF
  using CompoundSystemVariant =
      std::variant< legacy::resolved::CompoundSystem< SingleLevelBreitWigner >,
                    legacy::resolved::CompoundSystem< MultiLevelBreitWigner >,
                    CompoundSystem< ReichMoore, ShiftFactor >,
                    CompoundSystem< ReichMoore, Constant >,
                    //CompoundSystem< GeneralRMatrix, ShiftFactor >,
                    //CompoundSystem< GeneralRMatrix, Constant >,
                    legacy::unresolved::CompoundSystem >;
  #include "resonanceReconstruction/rmatrix/Reconstructor.hpp"
  #include "resonanceReconstruction/rmatrix/src/fromENDF.hpp"

} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
