#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction/rmatrix/CompoundSystem.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using Particle = rmatrix::Particle;
using ParticlePair = rmatrix::ParticlePair;
using ParticleID = rmatrix::ParticleID;
using ParticlePairID = rmatrix::ParticlePairID;
using Neutron = rmatrix::Neutron;
using Photon = rmatrix::Photon;
using Fission = rmatrix::Fission;
using ChargedParticle = rmatrix::ChargedParticle;
template < typename Type > using Channel = rmatrix::Channel< Type >;
using ParticleChannel = rmatrix::ParticleChannel;
using ParticleChannelData = rmatrix::ParticleChannelData;
using ChannelID = rmatrix::ChannelID;
using ChannelQuantumNumbers = rmatrix::ChannelQuantumNumbers;
using ChannelRadii = rmatrix::ChannelRadii;
using BoundaryCondition = rmatrix::BoundaryCondition;
using Resonance = rmatrix::Resonance;
using ResonanceTable = rmatrix::ResonanceTable;
template < typename Formalism, typename Option > using SpinGroup = rmatrix::SpinGroup< Formalism, Option >;
template < typename Formalism, typename Option > using CompoundSystem = rmatrix::CompoundSystem< Formalism, Option >;
using ReactionID = rmatrix::ReactionID;
using ReactionChannelID = rmatrix::ReactionChannelID;
using ShiftFactor = rmatrix::ShiftFactor;
using Constant = rmatrix::Constant;
using ReichMoore = rmatrix::ReichMoore;

constexpr AtomicMass neutronMass = 1.008664 * daltons;
constexpr ElectricalCharge elementary = dimwits::constant::elementaryCharge;

#include "resonanceReconstruction/rmatrix/CompoundSystem/test/CompoundSystem.test.hpp"
#include "resonanceReconstruction/rmatrix/CompoundSystem/test/evaluate.test.hpp"
#include "resonanceReconstruction/rmatrix/CompoundSystem/test/evaluateTMatrix.test.hpp"
