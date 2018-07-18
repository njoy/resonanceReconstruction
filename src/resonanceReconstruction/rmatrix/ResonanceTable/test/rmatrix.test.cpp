#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using ResonanceTable = rmatrix::ResonanceTable;
using Resonance = rmatrix::Resonance;
using ChannelID = rmatrix::ChannelID;

SCENARIO( "Resonance.rmatrix( energy )" ) {

  GIVEN( "Valid Resonances" ) {

    // three resonance table with two channels
    ResonanceTable table( { "1", "2" },
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
                                       7.451500e-1 * rootElectronVolt ) } );

    THEN( "an R-matrix contribution can be calculated" ) {

      auto rmatrix = table.rmatrix( 1e+4 * electronVolt );

      REQUIRE( rmatrix.rows() == 2 );
      REQUIRE( rmatrix.cols() == 2 );

      REQUIRE( 1.876492E+01 == Approx( rmatrix( 0, 0 ).real() ) );
      REQUIRE( 5.993566E-05 == Approx( rmatrix( 0, 0 ).imag() ) );
      REQUIRE( 4.080050E-03 == Approx( rmatrix( 0, 1 ).real() ) );
      REQUIRE( 1.313122E-08 == Approx( rmatrix( 0, 1 ).imag() ) );
      REQUIRE( 4.080050E-03 == Approx( rmatrix( 1, 0 ).real() ) );
      REQUIRE( 1.313122E-08 == Approx( rmatrix( 1, 0 ).imag() ) );
      REQUIRE( 9.274126E-07 == Approx( rmatrix( 1, 1 ).real() ) );
      REQUIRE( 9.274126E-07 == Approx( rmatrix( 1, 1 ).imag() ) );
    }
  } // GIVEN
} // SCENARIO


