/**
 *  @class
 *  @brief The l,S,J,pi quantum numbers of a reaction channel
 *
 *  The ChannelQuantumNumbers class contains the quantum numbers associated to 
 *  a given reaction channel. Only channels that have the same J,pi contribute
 *  to the cross section of a given reaction.
 */
class ChannelQuantumNumbers {

  /* fields */
  OrbitalAngularMomentum l_;
  Spin s_;
  TotalAngularMomentum J_;
  Parity parity_;

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/ChannelQuantumNumbers/src/ctor.hpp"

  /**
   *  @brief Return the orbital angular momentum l of the channel
   */
  const auto& orbitalAngularMomentum() const { return this->l_; }

  /**
   *  @brief Return the channel spin
   */
  const auto& spin() const { return this->s_; }

  /**
   *  @brief Return the total angular momentum J of the channel
   */
  const auto& totalAngularMomentum() const { return this->J_; }

  /**
   *  @brief Return the parity
   */
  const auto& parity() const { return this->parity_; }
};
