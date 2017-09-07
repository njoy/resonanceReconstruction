#include <vector>

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace dimwits;
using namespace njoy::resonanceReconstruction;

SCENARIO("radius"){
  SECTION("double overload"){
    const double AP = 10.0;
    
    auto reference = ranges::view::repeat_n( AP * rootBarns, 10 );

    auto trial = ranges::view::linear_distribute( -10, 30, 40 )
      | ranges::view::transform( []( double d ){ return std::pow( 2.0, d ); } )
      | ranges::view::transform( []( double d ){ return d * electronVolts ; } )
      | ranges::view::transform
        ( [ ap = radius( AP ) ]( auto&& e ){ return ap( e ); } );
    
    RANGES_FOR( const auto pair, ranges::view::zip( trial, reference ) ){
      const auto trial = std::get<0>(pair);
      const auto reference = std::get<1>(pair);
      REQUIRE( trial == reference );
    }
  }
  
  SECTION("TAB1 overload"){
    auto tab1 = []{
      auto regionPairs = std::make_tuple( std::vector< long >{ 3, 6 },
                                          std::vector< long >{ 1, 2 } );
      auto orderedPairs =
        std::make_tuple( std::vector< double >{ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 },
                         std::vector< double >{ 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 } );
  
      return njoy::ENDFtk::UnivariateTabulation( 1.0, 2.0, 3ul, 4ul,
                                                 std::move( regionPairs ),
                                                 std::move( orderedPairs ) );
    }();
    
    auto trial = ranges::view::linear_distribute( 1.0, 6.0, 11 )
      | ranges::view::transform( []( double e ){ return e * electronVolts; } )
      | ranges::view::transform
        ( [ ape = radius( tab1 ) ]( auto&& e ){ return ape( e ); } );
    
    std::vector< Quantity<RootBarn> > reference = { 3.0 * rootBarn,
                                                    3.0 * rootBarn,
                                                    4.0 * rootBarn,
                                                    4.0 * rootBarn,
                                                    5.0 * rootBarn,
                                                    5.5 * rootBarn,
                                                    6.0 * rootBarn,
                                                    6.5 * rootBarn,
                                                    7.0 * rootBarn,
                                                    7.5 * rootBarn,
                                                    8.0 * rootBarn };

    RANGES_FOR( const auto pair, ranges::view::zip( trial, reference ) ){
      const auto trial = std::get<0>(pair);
      const auto reference = std::get<1>(pair);
      REQUIRE( trial.value == Approx(reference.value) );
    }
  }
}