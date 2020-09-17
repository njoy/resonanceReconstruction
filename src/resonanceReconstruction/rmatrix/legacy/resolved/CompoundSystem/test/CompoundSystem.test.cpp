#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using Particle = rmatrix::Particle;
using ParticlePair = rmatrix::ParticlePair;
using ParticleID = rmatrix::ParticleID;
using ParticlePairID = rmatrix::ParticlePairID;
using Neutron = rmatrix::Neutron;
template < typename Type > using Channel = rmatrix::Channel< Type >;
using ChannelQuantumNumbers = rmatrix::ChannelQuantumNumbers;
using ChannelRadii = rmatrix::ChannelRadii;
using BoundaryCondition = rmatrix::BoundaryCondition;
using Resonance = rmatrix::legacy::resolved::Resonance;
using ResonanceTable = rmatrix::legacy::resolved::ResonanceTable;
template < typename Formalism >
using SpinGroup = rmatrix::legacy::resolved::SpinGroup< Formalism >;
template < typename Formalism >
using CompoundSystem = rmatrix::legacy::resolved::CompoundSystem< Formalism >;
using SingleLevelBreitWigner = rmatrix::legacy::resolved::SingleLevelBreitWigner;
using ReactionID = rmatrix::ReactionID;
using ReactionChannelID = rmatrix::ReactionChannelID;

constexpr AtomicMass neutronMass = 1.008664 * daltons;
constexpr ElectricalCharge elementary = dimwits::constant::elementaryCharge;
constexpr double e = 1.6021766208e-19;

#include "resonanceReconstruction/rmatrix/legacy/resolved/CompoundSystem/test/CompoundSystem.test.hpp"
#include "resonanceReconstruction/rmatrix/legacy/resolved/CompoundSystem/test/evaluate.test.hpp"
