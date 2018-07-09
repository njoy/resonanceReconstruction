namespace rmatrix {

  // quantum numbers
  using OrbitalAngularMomentum = unsigned int;
  using Spin = double;
  using TotalAngularMomentum = double;
  using Parity = short;

  // auxiliary functions for the quantum numbers
  #include "resonanceReconstruction/rmatrix/src/possibleChannelSpinValues.hpp"
  #include "resonanceReconstruction/rmatrix/src/possibleChannelTotalAngularMomentumValues.hpp"

  // R-Matrix boundary condition used in calculating the shift
  using BoundaryCondition = double;

  // identifiers
  using ReactionID = std::string;
  using ChannelID = std::string;

  #include "resonanceReconstruction/rmatrix/Particle.hpp"
  #include "resonanceReconstruction/rmatrix/ParticlePair.hpp"
  #include "resonanceReconstruction/rmatrix/ChannelQuantumNumbers.hpp"
  #include "resonanceReconstruction/rmatrix/Channel.hpp"

}
