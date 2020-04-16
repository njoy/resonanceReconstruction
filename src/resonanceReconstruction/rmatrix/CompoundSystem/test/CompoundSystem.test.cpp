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
using ReichMoore = rmatrix::ReichMoore;
template < typename Formalism, typename Option > using SpinGroup = rmatrix::SpinGroup< Formalism, Option >;
template < typename Formalism, typename Option > using CompoundSystem = rmatrix::CompoundSystem< Formalism, Option >;

constexpr AtomicMass neutronMass = 1.008664 * daltons;

SCENARIO( "CompoundSystem" ) {

  GIVEN( "valid data for a CompoundSystem" ) {

    // particles
    Particle photon( "g", 0.0 * daltons, 0.0 * coulombs, 1., +1);
    Particle neutron( "n", neutronMass, 0.0 * coulombs, 0.5, +1);
    Particle fe55( "Fe55_e0", 5.446635e+1 * neutronMass,
                              26.0 * coulombs, 0.0, +1);
    Particle fe54( "Fe54_e0", 5.347624e+1 * neutronMass,
                              26.0 * coulombs, 0.0, +1);

    // particle pairs
    ParticlePair out( photon, fe55, 0.0 * electronVolt );
    ParticlePair in( neutron, fe54, 0.0 * electronVolt );

    // channels
    Channel< Photon > capture1( out, { 0, 0.0, 0.5, +1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic1( in, { 0, 0.5, 0.5, +1 },
                                 { 5.437300e-1 * rootBarn,
                                   5.437300e-1 * rootBarn } );
    Channel< Photon > capture2( out, { 0, 0.0, 0.5, -1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic2( in, { 1, 0.5, 0.5, -1 },
                                 { 5.437300e-1 * rootBarn,
                                   5.437300e-1 * rootBarn } );
    Channel< Photon > capture3( out, { 0, 0.0, 1.5, -1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic3( in, { 1, 0.5, 1.5, -1 },
                                 { 5.437300e-1 * rootBarn,
                                   5.437300e-1 * rootBarn } );
    Channel< Photon > capture4( out, { 0, 0.0, 1.5, +1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic4( in, { 2, 0.5, 1.5, +1 },
                                 { 5.437300e-1 * rootBarn,
                                   5.437300e-1 * rootBarn } );
    Channel< Photon > capture5( out, { 0, 0.0, 2.5, +1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic5( in, { 2, 0.5, 2.5, +1 },
                                 { 5.437300e-1 * rootBarn,
                                   5.437300e-1 * rootBarn } );

    // conversion from Gamma to gamma
    auto eGamma = [&] ( double width, const Energy& energy,
                        const Channel< Neutron >& elastic ) -> ReducedWidth {
      return std::sqrt( width / 2. / elastic.penetrability( energy ) ) *
             rootElectronVolt;
    };
    auto cGamma = [&] ( double width ) -> ReducedWidth {
      return std::sqrt( width / 2. ) * rootElectronVolt;
    };

    // single resonance tables
    ResonanceTable table1(
      { elastic1.channelID() },
      { Resonance( 7.788000e+3 * electronVolt,
                   { eGamma( 1.187354e+3, 7.788000e+3 * electronVolt,
                             elastic1 ) },
                   cGamma( 1.455000e+0 ) ) } );
    ResonanceTable table2(
      { elastic2.channelID() },
      { Resonance( 5.152000e+4 * electronVolt,
                   { eGamma( 1.600200e+1, 5.152000e+4 * electronVolt,
                             elastic2 ) },
                   cGamma( 3.600000e-1 ) ) } );
    ResonanceTable table3(
      { elastic3.channelID() },
      { Resonance( 3.099000e+3 * electronVolt,
                   { eGamma( 1.400000e-3, 3.099000e+3 * electronVolt,
                             elastic3 ) },
                   cGamma( 5.900000e-1 ) ) } );
    ResonanceTable table4(
      { elastic4.channelID() },
      { Resonance( 9.480000e+3 * electronVolt,
                   { eGamma( 1.200000e+0, 9.480000e+3 * electronVolt,
                             elastic4 ) },
                   cGamma( 2.700000e-1 ) ) } );
    ResonanceTable table5(
      { elastic5.channelID() },
      { Resonance( 1.264000e+5 * electronVolt,
                   { eGamma( 2.900000e+0, 1.264000e+5 * electronVolt,
                             elastic5 ) },
                   cGamma( 1.100000e+0 ) ) } );

    SpinGroup< ReichMoore, ShiftFactor >
        group1( in, { elastic1 }, std::move( table1 ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group2( in, { elastic2 }, std::move( table2 ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group3( in, { elastic3 }, std::move( table3 ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group4( in, { elastic4 }, std::move( table4 ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group5( in, { elastic5 }, std::move( table5 ) );

    THEN( "a SpinGroup can be constructed" ) {
      CompoundSystem< ReichMoore, ShiftFactor >
          system( { group1, group2, group3, group4, group5 } );

      REQUIRE( 5 == system.spinGroups().size() );

      // group 1 - Jpi = 0.5+
      auto group = system.spinGroups()[0];
      REQUIRE( 1 == group.incidentChannels().size() );
      REQUIRE( 1 == group.channels().size() );
      auto channel =
          std::experimental::get< Channel< Neutron > >( group.channels()[0] );
      REQUIRE( "n,Fe54_e0{0,1/2,1/2+}" == channel.channelID() );

      // group 2 - Jpi = 0.5-
      group = system.spinGroups()[1];
      REQUIRE( 1 == group.incidentChannels().size() );
      REQUIRE( 1 == group.channels().size() );
      channel =
          std::experimental::get< Channel< Neutron > >( group.channels()[0] );
      REQUIRE( "n,Fe54_e0{1,1/2,1/2-}" == channel.channelID() );

      // group 3 - Jpi = 1.5-
      group = system.spinGroups()[2];
      REQUIRE( 1 == group.incidentChannels().size() );
      REQUIRE( 1 == group.channels().size() );
      channel =
          std::experimental::get< Channel< Neutron > >( group.channels()[0] );
      REQUIRE( "n,Fe54_e0{1,1/2,3/2-}" == channel.channelID() );

      // group 4 - Jpi = 1.5+
      group = system.spinGroups()[3];
      REQUIRE( 1 == group.incidentChannels().size() );
      REQUIRE( 1 == group.channels().size() );
      channel =
          std::experimental::get< Channel< Neutron > >( group.channels()[0] );
      REQUIRE( "n,Fe54_e0{2,1/2,3/2+}" == channel.channelID() );

      // group 5 - Jpi = 2.5+
      group = system.spinGroups()[4];
      REQUIRE( 1 == group.incidentChannels().size() );
      REQUIRE( 1 == group.channels().size() );
      channel =
          std::experimental::get< Channel< Neutron > >( group.channels()[0] );
      REQUIRE( "n,Fe54_e0{2,1/2,5/2+}" == channel.channelID() );
    } // THEN
  } // GIVEN

  GIVEN( "data for a CompoundSystem with errors" ) {

    // particles
    Particle photon( "g", 0.0 * daltons, 0.0 * coulombs, 1., +1);
    Particle neutron( "n", neutronMass, 0.0 * coulombs, 0.5, +1);
    Particle fe55( "Fe55_e0", 5.446635e+1 * neutronMass,
                              26.0 * coulombs, 0.0, +1);
    Particle fe54( "Fe54_e0", 5.347624e+1 * neutronMass,
                              26.0 * coulombs, 0.0, +1);

    // particle pairs
    ParticlePair out( photon, fe55, 0.0 * electronVolt );
    ParticlePair in( neutron, fe54, 0.0 * electronVolt );

    // channels
    Channel< Photon > capture1( out, { 0, 0.0, 0.5, +1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic1( in, { 0, 0.5, 0.5, +1 },
                                 { 5.437300e-1 * rootBarn,
                                   5.437300e-1 * rootBarn } );
    Channel< Photon > capture2( out, { 0, 0.0, 0.5, -1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic2( in, { 1, 0.5, 0.5, -1 },
                                 { 5.437300e-1 * rootBarn,
                                   5.437300e-1 * rootBarn } );
    Channel< Photon > capture3( out, { 0, 0.0, 1.5, -1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic3( in, { 1, 0.5, 1.5, -1 },
                                 { 5.437300e-1 * rootBarn,
                                   5.437300e-1 * rootBarn } );
    Channel< Photon > capture4( out, { 0, 0.0, 1.5, +1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic4( in, { 2, 0.5, 1.5, +1 },
                                 { 5.437300e-1 * rootBarn,
                                   5.437300e-1 * rootBarn } );

    // conversion from Gamma to gamma
    auto eGamma = [&] ( double width, const Energy& energy,
                        const Channel< Neutron >& elastic ) -> ReducedWidth {
      return std::sqrt( width / 2. / elastic.penetrability( energy ) ) *
             rootElectronVolt;
    };
    auto cGamma = [&] ( double width ) -> ReducedWidth {
      return std::sqrt( width / 2. ) * rootElectronVolt;
    };

    // multiple resonance table
    ResonanceTable table1(
      { elastic1.channelID() },
      { Resonance( 7.788000e+3 * electronVolt,
                   { eGamma( 1.187354e+3, 7.788000e+3 * electronVolt,
                             elastic1 ) },
                   cGamma( 1.455000e+0 ) ) } );
    ResonanceTable table2(
      { elastic2.channelID() },
      { Resonance( 5.152000e+4 * electronVolt,
                   { eGamma( 1.600200e+1, 5.152000e+4 * electronVolt,
                             elastic2 ) },
                   cGamma( 3.600000e-1 ) ) } );
    ResonanceTable table3(
      { elastic3.channelID() },
      { Resonance( 3.099000e+3 * electronVolt,
                   { eGamma( 1.400000e-3, 3.099000e+3 * electronVolt,
                             elastic3 ) },
                   cGamma( 5.900000e-1 ) ) } );
    ResonanceTable table4(
      { elastic4.channelID() },
      { Resonance( 9.480000e+3 * electronVolt,
                   { eGamma( 1.200000e+0, 9.480000e+3 * electronVolt,
                             elastic4 ) },
                   cGamma( 2.700000e-1 ) ) } );

    SpinGroup< ReichMoore, ShiftFactor >
        group1( in, { elastic1 }, std::move( table1 ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group2( in, { elastic2 }, std::move( table2 ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group3( in, { elastic3 }, std::move( table3 ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group4( in, { elastic4 }, std::move( table4 ) );

    // no need to test with other template parameter: construction of
    // CompoundSystem is independent of the template parameters
    using ReichMooreShiftFactorCompoundSystem =
        CompoundSystem< ReichMoore, ShiftFactor >;

    THEN( "an exception is thrown at construction when there are no spin "
          "groups" ) {

      REQUIRE_THROWS( ReichMooreShiftFactorCompoundSystem( {} ) );
    } // THEN

    THEN( "an exception is thrown at construction when the spin groups are not "
          "unique" ) {

      REQUIRE_THROWS( ReichMooreShiftFactorCompoundSystem
                          ( { group1, group2, group3, group4, group1 } ) );
    } // THEN
  } // GIVEN
} // SCENARIO


