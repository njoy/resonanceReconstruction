#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction/rmatrix/legacy/unresolved/CompoundSystem.hpp"

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
using Resonance = rmatrix::legacy::unresolved::Resonance;
using ResonanceTable = rmatrix::legacy::unresolved::ResonanceTable;
using SpinGroup = rmatrix::legacy::unresolved::SpinGroup;
using CompoundSystem = rmatrix::legacy::unresolved::CompoundSystem;
using ReactionID = rmatrix::ReactionID;
template < typename Key, typename Value > using Map = rmatrix::Map< Key, Value >;

constexpr AtomicMass neutronMass = 1.008664 * daltons;
constexpr ElectricalCharge elementary = dimwits::constant::elementaryCharge;
constexpr double e = 1.6021766208e-19;

#include "resonanceReconstruction/rmatrix/legacy/unresolved/CompoundSystem/test/CompoundSystem.test.hpp"
#include "resonanceReconstruction/rmatrix/legacy/unresolved/CompoundSystem/test/evaluate.test.hpp"
