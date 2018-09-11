#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using ResonanceTable = rmatrix::ResonanceTable;
using Resonance = rmatrix::Resonance;
using ChannelID = rmatrix::ChannelID;

SCENARIO( "ResonanceTable" ) {

  GIVEN( "valid data for a ResonanceTable" ) {

    // single resonance data
    std::vector< ChannelID > channels = { "1", "2" };
    std::vector< Resonance > resonances =
      { Resonance( 6.823616e+4 * electronVolt,
                   { 2.179040e+2 * rootElectronVolt,
                     1.000000e-5 * rootElectronVolt }, 
                   3.933600e-1 * rootElectronVolt ),
        Resonance( 1.150980e+5 * electronVolt,
                   { 4.307780e+0 * rootElectronVolt,
                     0.0 * rootElectronVolt }, 
                   7.390000e-1 * rootElectronVolt ),
        Resonance( 1.825230e+5 * electronVolt,
                   { 1.759740e+3 * rootElectronVolt,
                     4.000000e-1 * rootElectronVolt }, 
                   7.451500e-1 * rootElectronVolt ) };

    THEN( "a ResonanceTable can be constructed" ) {
      ResonanceTable table( std::move( channels ), std::move( resonances ) );

      REQUIRE( 2 == table.numberChannels() );
      REQUIRE( 2 == table.channels().size() );
      REQUIRE( "1" == table.channels()[0] );
      REQUIRE( "2" == table.channels()[1] );

      REQUIRE( 3 == table.numberResonances() );
      REQUIRE( 3 == table.resonances().size() );
      REQUIRE( 3 == table.energies().size() );
      REQUIRE( 6.823616e+4 == Approx( table.energies()[0].value ) );
      REQUIRE( 1.150980e+5 == Approx( table.energies()[1].value ) );
      REQUIRE( 1.825230e+5 == Approx( table.energies()[2].value ) );

      auto resonance = table.resonances()[0];
      REQUIRE( 6.823616e+4 == Approx( resonance.energy().value ) );
      REQUIRE( 3.933600e-1 == Approx( resonance.eliminatedWidth().value ) );
      REQUIRE( 2 == resonance.widths().size() );
      REQUIRE( 2.179040e+2 == Approx( resonance.widths()[0].value ) );
      REQUIRE( 1.000000e-5 == Approx( resonance.widths()[1].value ) );

      resonance = table.resonances()[1];
      REQUIRE( 1.150980e+5 == Approx( resonance.energy().value ) );
      REQUIRE( 7.390000e-1 == Approx( resonance.eliminatedWidth().value ) );
      REQUIRE( 2 == resonance.widths().size() );
      REQUIRE( 4.307780e+0 == Approx( resonance.widths()[0].value ) );
      REQUIRE( 0.0 == Approx( resonance.widths()[1].value ) );

      resonance = table.resonances()[2];
      REQUIRE( 1.825230e+5 == Approx( resonance.energy().value ) );
      REQUIRE( 7.451500e-1 == Approx( resonance.eliminatedWidth().value ) );
      REQUIRE( 2 == resonance.widths().size() );
      REQUIRE( 1.759740e+3 == Approx( resonance.widths()[0].value ) );
      REQUIRE( 4.000000e-1 == Approx( resonance.widths()[1].value ) );
    }
  } // GIVEN

  GIVEN( "data for a ResonanceTable containing errors" ) {

    // wrong number of widths in a Resonance
    std::vector< ChannelID > channels = { "1", "2" };
    std::vector< Resonance > resonances =
      { Resonance( 6.823616e+4 * electronVolt,
                   { 2.179040e+2 * rootElectronVolt,
                     1.000000e-5 * rootElectronVolt },
                   3.933600e-1 * rootElectronVolt ),
        Resonance( 1.150980e+5 * electronVolt,
                   { 4.307780e+0 * rootElectronVolt },
                   7.390000e-1 * rootElectronVolt ),
        Resonance( 1.825230e+5 * electronVolt,
                   { 1.759740e+3 * rootElectronVolt,
                     4.000000e-1 * rootElectronVolt },
                   7.451500e-1 * rootElectronVolt ) };

    THEN( "an exception is thrown at construction when a resonance does not "
          "contain the right number of widths" ) {

      REQUIRE_THROWS( ResonanceTable( std::move( channels ),
                                      std::move( resonances ) ) );
    }
  } // GIVEN
} // SCENARIO


