#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
template < typename Quantity >
using Data = rmatrix::legacy::Data< Quantity >;

SCENARIO( "Data" ) {

  GIVEN( "valid data for a Data object" ) {

    THEN( "a Data object can be constructed" ) {

      Data< double > doubles( 1., 2., 3., 4. );

      CHECK( 1.0 == doubles.elastic );
      CHECK( 2.0 == doubles.capture );
      CHECK( 3.0 == doubles.fission );
      CHECK( 4.0 == doubles.competition );

      CHECK( true == doubles.hasElastic() );
      CHECK( true == doubles.hasCapture() );
      CHECK( true == doubles.hasFission() );
      CHECK( true == doubles.hasCompetition() );

      Data< double > zeros( 0., 0., 0., 0. );

      CHECK( 0.0 == zeros.elastic );
      CHECK( 0.0 == zeros.capture );
      CHECK( 0.0 == zeros.fission );
      CHECK( 0.0 == zeros.competition );

      CHECK( false == zeros.hasElastic() );
      CHECK( false == zeros.hasCapture() );
      CHECK( false == zeros.hasFission() );
      CHECK( false == zeros.hasCompetition() );

      Data< double > some( 1., 2., 0., 0. );

      CHECK( 1.0 == some.elastic );
      CHECK( 2.0 == some.capture );
      CHECK( 0.0 == some.fission );
      CHECK( 0.0 == some.competition );

      CHECK( true == some.hasElastic() );
      CHECK( true == some.hasCapture() );
      CHECK( false == some.hasFission() );
      CHECK( false == some.hasCompetition() );
    } // THEN
  } // GIVEN

  GIVEN( "valid Data objects" ) {

    Data< double > data1( 1., 2., 3., 4. );
    Data< double > data2( 4., 3., 2., 1. );

    THEN( "operator overloading can be used" ) {

      Data< double > result = data1 + data2;
      CHECK( 5.0 == result.elastic );
      CHECK( 5.0 == result.capture );
      CHECK( 5.0 == result.fission );
      CHECK( 5.0 == result.competition );

      result = data1 - data2;
      CHECK( -3.0 == result.elastic );
      CHECK( -1.0 == result.capture );
      CHECK( 1.0 == result.fission );
      CHECK( 3.0 == result.competition );

      result += data2;
      CHECK( 1.0 == result.elastic );
      CHECK( 2.0 == result.capture );
      CHECK( 3.0 == result.fission );
      CHECK( 4.0 == result.competition );

      result -= data1;
      CHECK( 0.0 == result.elastic );
      CHECK( 0.0 == result.capture );
      CHECK( 0.0 == result.fission );
      CHECK( 0.0 == result.competition );
    } // THEN
  } // GIVEN
} // SCENARIO
