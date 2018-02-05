#define CATCH_CONFIG_MAIN

#include <iterator>
#include <fstream>
#include <tuple>

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;
using namespace dimwits;

std::pair< njoy::ENDFtk::section::Type<2>, std::vector< double > >
resonances( const std::string& id );

static void escape( void* );

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

auto timingTestInclusive( const std::vector< double >& testData ){
  return [ &testData ]( auto&& xs ){
    auto energies = testData | ranges::view::stride( 4 );
    
    auto accumulation = 0.0*barns;
    RANGES_FOR( auto energy, energies ){
      accumulation += xs( energy*electronVolts ).elastic;
    }
    escape( &accumulation );
  };
}

auto timingTestExclusive( const std::vector< double >& testData ){
  return [ &testData ]( auto&& xs ){
    auto energies = testData | ranges::view::stride( 4 );
    
    auto start = std::chrono::high_resolution_clock::now();

    auto accumulation = 0.0*barns;
    RANGES_FOR( auto energy, energies ){
      accumulation += xs( energy*electronVolts ).elastic;
    }
    escape( &accumulation );

    auto finish = std::chrono::high_resolution_clock::now();
    auto microseconds =
      std::chrono::duration_cast<std::chrono::microseconds>(finish-start);
    njoy::Log::info( "Approximately {} microseconds passed while 'Apply'-ing",
                    microseconds.count() );
  };
}

SCENARIO( "Integration test" ){
  SECTION( "Cobalt-58" ){
    auto Co58 = resonances("Co-58");

    auto& section151 = std::get<0>( Co58 );
    auto& isotope = section151.isotopes.front();
    auto& energyRange = isotope.energyRanges().front();
    auto& mlbw = std::experimental::get< 2 >( energyRange );

    njoy::Log::info("\n Cobalt-58 "
                    "\n --------------- "
                    "\n LRU: {}"
                    "\n LRF: {}"
                    "\n NRO: {}"
                    "\n NAPS: {}\n",
                    mlbw.LRU(), mlbw.LRF(), mlbw.NRO(), mlbw.NAPS() );

    auto& testData = std::get<1>( Co58 );

    breitWigner::multilevel::Apply{}( mlbw, test( testData ) );
  }

  SECTION( "Tulium-168" ){
    auto Tm168 = resonances("Tm-168");

    auto& section151 = std::get<0>( Tm168 );
    auto& isotope = section151.isotopes.front();
    auto& energyRange = isotope.energyRanges().front();
    auto& mlbw = std::experimental::get< 2 >( energyRange );

    njoy::Log::info("\n Tulium-168 "
                    "\n --------------- "
                    "\n LRU: {}"
                    "\n LRF: {}"
                    "\n NRO: {}"
                    "\n NAPS: {}\n",
                    mlbw.LRU(), mlbw.LRF(), mlbw.NRO(), mlbw.NAPS() );

    auto& testData = std::get<1>( Tm168 );

    breitWigner::multilevel::Apply{}( mlbw, test( testData ) );
  }

  SECTION( "Neptunium-238" ){
    auto Np238 = resonances("Np-238");

    auto& section151 = std::get<0>( Np238 );
    auto& isotope = section151.isotopes.front();
    auto& energyRange = isotope.energyRanges().front();
    auto& mlbw = std::experimental::get< 2 >( energyRange );

    njoy::Log::info( "\n Neptunium-238 "
                     "\n --------------- "
                     "\n LRU: {}"
                     "\n LRF: {}"
                     "\n NRO: {}"
                     "\n NAPS: {}\n",
                     mlbw.LRU(), mlbw.LRF(), mlbw.NRO(), mlbw.NAPS() );

    auto& testData = std::get<1>( Np238 );

    breitWigner::multilevel::Apply{}( mlbw, test( testData ) );
  }
}

SCENARIO( "Timing test" ){
  SECTION( "Cobalt-58" ){
    auto Co58 = resonances("Co-58");

    auto& section151 = std::get<0>( Co58 );
    auto& isotope = section151.isotopes.front();
    auto& energyRange = isotope.energyRanges().front();
    auto& mlbw = std::experimental::get< 2 >( energyRange );

    njoy::Log::info("\n Cobalt-58 "
                    "\n --------------- ");

    auto& testData = std::get<1>( Co58 );

    SECTION( "\nInclusive timing" ){
      njoy::Log::info( "inclusive timing" );
      auto start = std::chrono::high_resolution_clock::now();

      breitWigner::multilevel::Apply{}( mlbw, timingTestInclusive( testData ) );

      auto finish = std::chrono::high_resolution_clock::now();
      auto microseconds =
        std::chrono::duration_cast<std::chrono::microseconds>(finish-start);
      njoy::Log::info( "Approximately {} microseconds passed while 'Apply'-ing",
                      microseconds.count() );
    }

    SECTION( "\nExclusive timing" ){
      njoy::Log::info( "exclusive timing" );
      breitWigner::multilevel::Apply{}( mlbw, timingTestExclusive( testData ) );
    }
  }

  SECTION( "Tulium-168" ){
    auto Tm168 = resonances("Tm-168");

    auto& section151 = std::get<0>( Tm168 );
    auto& isotope = section151.isotopes.front();
    auto& energyRange = isotope.energyRanges().front();
    auto& mlbw = std::experimental::get< 2 >( energyRange );

    njoy::Log::info("\n Tulium-168 "
                    "\n --------------- ");

    auto& testData = std::get<1>( Tm168 );

    SECTION( "\nInclusive timing" ){
      njoy::Log::info( "inclusive timing" );
      auto start = std::chrono::high_resolution_clock::now();

      breitWigner::multilevel::Apply{}( mlbw, timingTestInclusive( testData ) );

      auto finish = std::chrono::high_resolution_clock::now();
      auto microseconds =
        std::chrono::duration_cast<std::chrono::microseconds>(finish-start);
      njoy::Log::info( "Approximately {} microseconds passed while 'Apply'-ing",
                      microseconds.count() );
    }

    SECTION( "\nExclusive timing" ){
      njoy::Log::info( "exclusive timing" );
      breitWigner::multilevel::Apply{}( mlbw, timingTestExclusive( testData ) );
    }
  }

  SECTION( "Neptunium-238" ){
    auto Np238 = resonances("Np-238");

    auto& section151 = std::get<0>( Np238 );
    auto& isotope = section151.isotopes.front();
    auto& energyRange = isotope.energyRanges().front();
    auto& mlbw = std::experimental::get< 2 >( energyRange );

    njoy::Log::info( "\n Neptunium-238 "
                     "\n --------------- ");

    auto& testData = std::get<1>( Np238 );

    SECTION( "\nInclusive timing" ){
      njoy::Log::info( "inclusive timing" );
      auto start = std::chrono::high_resolution_clock::now();

      breitWigner::multilevel::Apply{}( mlbw, timingTestInclusive( testData ) );

      auto finish = std::chrono::high_resolution_clock::now();
      auto microseconds =
        std::chrono::duration_cast<std::chrono::microseconds>(finish-start);
      njoy::Log::info( "Approximately {} microseconds passed while 'Apply'-ing",
                      microseconds.count() );
    }

    SECTION( "\nExclusive timing" ){
      njoy::Log::info( "exclusive timing" );
      breitWigner::multilevel::Apply{}( mlbw, timingTestExclusive( testData ) );
    }
  }
}

std::pair< njoy::ENDFtk::section::Type<2>, std::vector< double > >
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

    auto begin = endfFile.begin();
    auto end = endfFile.end();

    njoy::ENDFtk::syntaxTree::Tape< std::string::iterator > tape( begin, end );

    auto& material = *( tape.begin() );

    auto MAT = material.MAT();
    long lineNumber = 1;
    return material
           .fileNumber(2)
           .sectionNumber(151).parse<2>();
  };

  return std::make_pair( section151(), testData() );
}

// Force the compiler to not-eliminate some piece of memory
// Thus forcing the compiler to not throw away it's construction
static void escape( void* p ){
  asm volatile( "" : : "g"(p) : "memory" );
}
