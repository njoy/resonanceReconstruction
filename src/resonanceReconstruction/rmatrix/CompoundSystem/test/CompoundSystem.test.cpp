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
using Sammy = rmatrix::Sammy;
template < typename Option > using SpinGroup = rmatrix::SpinGroup< Option >;
template < typename Option > using CompoundSystem = rmatrix::CompoundSystem< Option >;

SCENARIO( "CompoundSystem" ) {

  GIVEN( "valid data for a CompoundSystem" ) {

  } // GIVEN
} // SCENARIO


