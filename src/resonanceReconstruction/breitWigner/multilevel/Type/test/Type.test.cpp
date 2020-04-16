#define CATCH_CONFIG_MAIN

#include <iterator>
#include <fstream>
#include <tuple>

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;
using namespace dimwits;

std::pair< njoy::ENDFtk::section::Type< 2, 151 >, std::vector< double > >
resonances( const std::string& id );

auto test( const std::vector< double >& testData ){
  return [&testData]( auto&& xs ){
    auto tuples = testData | ranges::view::chunk(4);
    RANGES_FOR( auto tuple, tuples ){
      auto energy = tuple[0] * electronVolts;
      double referenceElastic = tuple[1];
      double referenceFission = tuple[2];
      double referenceCapture = tuple[3];

      auto trial = xs( energy );
      REQUIRE( referenceElastic == Approx( trial.elastic.value ) );
      REQUIRE( referenceFission == Approx( trial.fission.value ) );
      REQUIRE( referenceCapture == Approx( trial.capture.value ) );
    }
  };
}

SCENARIO( "Integration test" ){
  SECTION( "Cobalt-58" ){
    auto Co58 = resonances("Co-58");

    auto& section151 = std::get<0>( Co58 );
    auto& isotope = section151.isotopes().front();
    auto& resonanceRange = isotope.resonanceRanges().front();

    njoy::Log::info("\n Cobalt-58 "
                    "\n --------------- "
                    "\n LRU: {}"
                    "\n LRF: {}"
                    "\n NRO: {}"
                    "\n NAPS: {}\n",
                    resonanceRange.LRU(), resonanceRange.LRF(),
                    resonanceRange.NRO(), resonanceRange.NAPS() );

    auto& testData = std::get<1>( Co58 );

    breitWigner::multilevel::Apply{}( resonanceRange, test( testData ) );
  }

  /*
  SECTION( "Tulium-168" ){
    auto Tm168 = resonances("Tm-168");

    auto& section151 = std::get<0>( Tm168 );
    auto& isotope = section151.isotopes().front();
    auto& resonanceRange = isotope.resonanceRanges().front();

    njoy::Log::info("\n Tulium-168 "
                    "\n --------------- "
                    "\n LRU: {}"
                    "\n LRF: {}"
                    "\n NRO: {}"
                    "\n NAPS: {}\n",
                    resonanceRange.LRU(), resonanceRange.LRF(),
                    resonanceRange.NRO(), resonanceRange.NAPS() );

    auto& testData = std::get<1>( Tm168 );

    breitWigner::multilevel::Apply{}( resonanceRange, test( testData ) );
  }

  SECTION( "Neptunium-238" ){
    auto Np238 = resonances("Np-238");

    auto& section151 = std::get<0>( Np238 );
    auto& isotope = section151.isotopes().front();
    auto& resonanceRange = isotope.resonanceRanges().front();

    njoy::Log::info( "\n Neptunium-238 "
                     "\n --------------- "
                     "\n LRU: {}"
                     "\n LRF: {}"
                     "\n NRO: {}"
                     "\n NAPS: {}\n",
                     resonanceRange.LRU(), resonanceRange.LRF(),
                     resonanceRange.NRO(), resonanceRange.NAPS() );

    auto& testData = std::get<1>( Np238 );

    breitWigner::multilevel::Apply{}( resonanceRange, test( testData ) );
  }
  */
}

std::pair< njoy::ENDFtk::section::Type< 2, 151 >, std::vector< double > >
resonances( const std::string& id ){
  auto testData = [&]{
    std::vector< double > data;
    std::ifstream tupleFile( id + "-tuples.txt" );
    for ( std::istream_iterator< double > it( tupleFile );
          it != std::istream_iterator< double >();
          ++it ){
      data.push_back( *it );
    }
    return data;
  };

  auto section151 = [&]{

    auto endfFile = njoy::utility::slurpFileToMemory( id + ".endf" );

    njoy::ENDFtk::syntaxTree::Tape< std::string > tape( endfFile );

    auto& material = *( tape.begin() );

    auto MAT = material.MAT();
    long lineNumber = 1;
    return material
           .fileNumber(2)
           .sectionNumber(151).parse< 2, 151 >();
  };

  return std::make_pair( section151(), testData() );
}
