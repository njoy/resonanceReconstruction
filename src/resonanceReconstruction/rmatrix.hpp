// physical quantities
#include "resonanceReconstruction/Quantity.hpp"

// formalism options, boundary condition options and channel types
#include "resonanceReconstruction/rmatrix/Formalism.hpp"
#include "resonanceReconstruction/rmatrix/BoundaryOption.hpp"
#include "resonanceReconstruction/rmatrix/ChannelType.hpp"

// R-matrix components (independent of formalism)
#include "resonanceReconstruction/rmatrix/QuantumNumbers.hpp"
#include "resonanceReconstruction/rmatrix/ChannelQuantumNumbers.hpp"
#include "resonanceReconstruction/rmatrix/Particle.hpp"
#include "resonanceReconstruction/rmatrix/ParticlePair.hpp"
#include "resonanceReconstruction/rmatrix/ChannelID.hpp"
#include "resonanceReconstruction/rmatrix/ChannelRadiusTable.hpp"
#include "resonanceReconstruction/rmatrix/ChannelRadii.hpp"
#include "resonanceReconstruction/rmatrix/Channel.hpp"
#include "resonanceReconstruction/rmatrix/ParticleChannel.hpp"
#include "resonanceReconstruction/rmatrix/ParticleChannelData.hpp"
#include "resonanceReconstruction/rmatrix/ReactionChannelID.hpp"

// resonance information
#include "resonanceReconstruction/rmatrix/Resonance.hpp"
#include "resonanceReconstruction/rmatrix/ResonanceTable.hpp"

// option dependent calculators
#include "resonanceReconstruction/rmatrix/LMatrixCalculator.hpp"
#include "resonanceReconstruction/rmatrix/RLMatrixCalculator.hpp"

// spin group and compound system
#include "resonanceReconstruction/rmatrix/SpinGroup.hpp"
#include "resonanceReconstruction/rmatrix/CompoundSystem.hpp"

// legacy resonance reconstruction
#include "resonanceReconstruction/rmatrix/legacy.hpp"

// the reconstructor object anf the fromENDF function
#include "resonanceReconstruction/rmatrix/Reconstructor.hpp"
#include "resonanceReconstruction/rmatrix/fromENDF.hpp"

// wave function calculation
#include "resonanceReconstruction/rmatrix/calculatePenetrability.hpp"
#include "resonanceReconstruction/rmatrix/calculateShiftFactor.hpp"
#include "resonanceReconstruction/rmatrix/calculatePhaseShift.hpp"
#include "resonanceReconstruction/rmatrix/calculateCoulombPhaseShift.hpp"

// auxiliary functions for the quantum numbers
#include "resonanceReconstruction/rmatrix/possibleChannelSpinValues.hpp"
#include "resonanceReconstruction/rmatrix/possibleChannelTotalAngularMomentumValues.hpp"
