#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction/rmatrix/ChannelRadiusTable.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using RadiusTable = rmatrix::MultiRegionTable< Energy, ChannelRadius >;
using HistogramTable = rmatrix::HistogramTable< Energy, ChannelRadius >;
using LinLinTable = rmatrix::LinLinTable< Energy, ChannelRadius >;
using LinLogTable = rmatrix::LinLogTable< Energy, ChannelRadius >;
using LogLinTable = rmatrix::LogLinTable< Energy, ChannelRadius >;
using LogLogTable = rmatrix::LogLogTable< Energy, ChannelRadius >;
using TableVariant = rmatrix::TableVariant< Energy, ChannelRadius >;
using ChannelRadiusTable = rmatrix::ChannelRadiusTable;

SCENARIO( "ChannelRadiusTable" ) {

  GIVEN( "valid data for a ChannelRadiusTable" ) {

    std::vector< Energy > energies = { 1e-5 * electronVolt, 2e+7 * electronVolt };
    std::vector< ChannelRadius > radii = { .5 * rootBarn, 1.0 * rootBarn };
    RadiusTable singletable( std::vector< TableVariant >{
                   LinLinTable( std::move( energies ), std::move( radii ) ) } );

    std::vector< Energy > energies1 = { 1e-5 * electronVolt, 1.0 * electronVolt };
    std::vector< Energy > energies2 = { 1. * electronVolt, 2e+7 * electronVolt };
    std::vector< ChannelRadius > radii1 = { .5 * rootBarn, 1.0 * rootBarn };
    std::vector< ChannelRadius > radii2 = { 1. * rootBarn, 1. * rootBarn };
    RadiusTable multitable( std::vector< TableVariant >{
                   LinLinTable( std::move( energies1 ), std::move( radii1 ) ),
                   HistogramTable( std::move( energies2 ), std::move( radii2 ) ) } );

    THEN( "a ChannelRadiusTable can be constructed for a single region" ) {

      ChannelRadiusTable radius( std::move( singletable ) );

      CHECK( .5 == Approx( radius( 1e-5 * electronVolt ).value ) );
      CHECK( 5.000000000025E-01 == Approx( radius( 1e-4 * electronVolt ).value ) );
      CHECK( 5.000000000250E-01 == Approx( radius( 1e-3 * electronVolt ).value ) );
      CHECK( 5.000000002500E-01 == Approx( radius( 1e-2 * electronVolt ).value ) );
      CHECK( 5.000000025000E-01 == Approx( radius( 1e-1 * electronVolt ).value ) );
      CHECK( 5.000000250000E-01 == Approx( radius( 1. * electronVolt ).value ) );
      CHECK( 5.000002500000E-01 == Approx( radius( 1e+1 * electronVolt ).value ) );
      CHECK( 5.000025000000E-01 == Approx( radius( 1e+2 * electronVolt ).value ) );
      CHECK( 5.000250000000E-01 == Approx( radius( 1e+3 * electronVolt ).value ) );
      CHECK( 5.002500000000E-01 == Approx( radius( 1e+4 * electronVolt ).value ) );
      CHECK( 5.025000000000E-01 == Approx( radius( 1e+5 * electronVolt ).value ) );
      CHECK( 5.250000000000E-01 == Approx( radius( 1e+6 * electronVolt ).value ) );
      CHECK( 7.500000000000E-01 == Approx( radius( 1e+7 * electronVolt ).value ) );
      CHECK( 1. == Approx( radius( 2e+7 * electronVolt ).value ) );
    } // THEN

    THEN( "a ChannelRadiusTable can be constructed for multiple regions" ) {

      ChannelRadiusTable radius( std::move( multitable ) );

      CHECK( .5 == Approx( radius( 1e-5 * electronVolt ).value ) );
      CHECK( 5.00045000450E-01 == Approx( radius( 1e-4 * electronVolt ).value ) );
      CHECK( 5.00495004950E-01 == Approx( radius( 1e-3 * electronVolt ).value ) );
      CHECK( 5.04995049950E-01 == Approx( radius( 1e-2 * electronVolt ).value ) );
      CHECK( 5.49995499955E-01 == Approx( radius( 1e-1 * electronVolt ).value ) );
      CHECK( 1. == Approx( radius( 1. * electronVolt ).value ) );
      CHECK( 1. == Approx( radius( 1e+1 * electronVolt ).value ) );
      CHECK( 1. == Approx( radius( 1e+2 * electronVolt ).value ) );
      CHECK( 1. == Approx( radius( 1e+3 * electronVolt ).value ) );
      CHECK( 1. == Approx( radius( 1e+4 * electronVolt ).value ) );
      CHECK( 1. == Approx( radius( 1e+5 * electronVolt ).value ) );
      CHECK( 1. == Approx( radius( 1e+6 * electronVolt ).value ) );
      CHECK( 1. == Approx( radius( 1e+7 * electronVolt ).value ) );
      CHECK( 1. == Approx( radius( 2e+7 * electronVolt ).value ) );
    } // THEN
  } // GIVEN
} // SCENARIO
