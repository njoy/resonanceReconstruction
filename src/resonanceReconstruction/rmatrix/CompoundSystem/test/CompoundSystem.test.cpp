#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using Particle = rmatrix::Particle;
using ParticlePair = rmatrix::ParticlePair;
using Neutron = rmatrix::Neutron;
using Photon = rmatrix::Photon;
using Fission = rmatrix::Fission;
using ReactionID = rmatrix::ReactionID;
template < typename Type > using Channel = rmatrix::Channel< Type >;
using Resonance = rmatrix::Resonance;
using ResonanceTable = rmatrix::ResonanceTable;
using ShiftFactor = rmatrix::ShiftFactor;
using Constant = rmatrix::Constant;
template < typename Formalism, typename Option > using SpinGroup = rmatrix::SpinGroup< Formalism, Option >;
template < typename Formalism, typename Option > using CompoundSystem = rmatrix::CompoundSystem< Formalism, Option >;

SCENARIO( "CompoundSystem" ) {

  GIVEN( "valid data for a CompoundSystem" ) {

  } // GIVEN
} // SCENARIO


