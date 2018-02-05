#define CATCH_CONFIG_MAIN

#include <chrono>

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
    auto milliseconds =
      std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
    njoy::Log::info( "Approximately {} milliseconds passed while 'Apply'-ing",
                    milliseconds.count() );
  };
}

SCENARIO( "Timing test" ){
  SECTION( "Iron-56" ){
    const auto Fe56 = resonances("Fe-56");

    const auto& section151 = std::get<0>( Fe56 );
    const auto& isotope = section151.isotopes.front();
    const auto& energyRange = isotope.energyRanges().front();
    const auto& rm = std::experimental::get< 3 >( energyRange );

    njoy::Log::info("\n Iron-56"
                    "\n --------------"
                    "\n LRU: {}"
                    "\n LRF: {}"
                    "\n NRO: {}"
                    "\n NAPS: {}\n",
                    rm.LRU(), rm.LRF(), rm.NRO(), rm.NAPS() );

    auto& testData = std::get<1>( Fe56 );

    SECTION( "Inclusive" ){
      njoy::Log::info( "inclusive timing" );
      auto start = std::chrono::high_resolution_clock::now();

      reichMoore::Apply{}( rm, timingTestInclusive( testData ) );

      auto finish = std::chrono::high_resolution_clock::now();
      auto milliseconds =
        std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
      njoy::Log::info( "Approximately {} milliseconds passed while 'Apply'-ing",
                      milliseconds.count() );
    }
    
    SECTION( "Exclusive" ){
      njoy::Log::info( "exclusive timing" );
      reichMoore::Apply{}( rm, timingTestExclusive( testData ) );
    }
  }

  SECTION( "Uranium-235" ){
    const auto U235 = resonances("U-235");

    const auto& section151 = std::get<0>( U235 );
    const auto& isotope = section151.isotopes.front();
    const auto& energyRange = isotope.energyRanges().front();
    const auto& rm = std::experimental::get< 3 >( energyRange );

    njoy::Log::info("\n Uranium-235"
                    "\n --------------"
                    "\n LRU: {}"
                    "\n LRF: {}"
                    "\n NRO: {}"
                    "\n NAPS: {}\n",
                    rm.LRU(), rm.LRF(), rm.NRO(), rm.NAPS() );

    auto& testData = std::get<1>( U235 );

    SECTION( "Inclusive" ){
      njoy::Log::info( "inclusive timing" );
      auto start = std::chrono::high_resolution_clock::now();

      reichMoore::Apply{}( rm, timingTestInclusive( testData ) );

      auto finish = std::chrono::high_resolution_clock::now();
      auto milliseconds =
        std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
      njoy::Log::info( "Approximately {} milliseconds passed while 'Apply'-ing",
                      milliseconds.count() );
    }
    
    SECTION( "Exclusive" ){
      njoy::Log::info( "exclusive timing" );
      reichMoore::Apply{}( rm, timingTestExclusive( testData ) );
    }
  }

  SECTION( "Uranium-238" ){
    const auto U238 = resonances("U-238");

    const auto& section151 = std::get<0>( U238 );
    const auto& isotope = section151.isotopes.front();
    const auto& energyRange = isotope.energyRanges().front();
    const auto& rm = std::experimental::get< 3 >( energyRange );

    njoy::Log::info("\n Uranium-238"
                    "\n --------------"
                    "\n LRU: {}"
                    "\n LRF: {}"
                    "\n NRO: {}"
                    "\n NAPS: {}\n",
                    rm.LRU(), rm.LRF(), rm.NRO(), rm.NAPS() );

    auto& testData = std::get<1>( U238 );

    SECTION( "Inclusive" ){
      njoy::Log::info( "inclusive timing" );
      auto start = std::chrono::high_resolution_clock::now();

      reichMoore::Apply{}( rm, timingTestInclusive( testData ) );

      auto finish = std::chrono::high_resolution_clock::now();
      auto milliseconds =
        std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
      njoy::Log::info( "Approximately {} milliseconds passed while 'Apply'-ing",
                      milliseconds.count() );
    }
    
    SECTION( "Exclusive" ){
      njoy::Log::info( "exclusive timing" );
      reichMoore::Apply{}( rm, timingTestExclusive( testData ) );
    }
  }
}

SCENARIO( "Integration test" ){
  SECTION( "Iron-56" ){
    const auto Fe56 = resonances("Fe-56");

    const auto& section151 = std::get<0>( Fe56 );
    const auto& isotope = section151.isotopes.front();
    const auto& energyRange = isotope.energyRanges().front();
    const auto& rm = std::experimental::get< 3 >( energyRange );

    njoy::Log::info("\n Iron-56"
                    "\n --------------"
                    "\n LRU: {}"
                    "\n LRF: {}"
                    "\n NRO: {}"
                    "\n NAPS: {}\n",
                    rm.LRU(), rm.LRF(), rm.NRO(), rm.NAPS() );

    auto& testData = std::get<1>( Fe56 );
    reichMoore::Apply{}( rm, test( testData ) );
  }

  SECTION( "Uranium-235" ){
    const auto U235 = resonances("U-235");

    const auto& section151 = std::get<0>( U235 );
    const auto& isotope = section151.isotopes.front();
    const auto& energyRange = isotope.energyRanges().front();
    const auto& rm = std::experimental::get< 3 >( energyRange );

    njoy::Log::info("\n Uranium-235"
                    "\n --------------"
                    "\n LRU: {}"
                    "\n LRF: {}"
                    "\n NRO: {}"
                    "\n NAPS: {}\n",
                    rm.LRU(), rm.LRF(), rm.NRO(), rm.NAPS() );

    auto& testData = std::get<1>( U235 );
    reichMoore::Apply{}( rm, test( testData ) );
  }

  SECTION( "Uranium-238" ){
    const auto U238 = resonances("U-238");

    const auto& section151 = std::get<0>( U238 );
    const auto& isotope = section151.isotopes.front();
    const auto& energyRange = isotope.energyRanges().front();
    const auto& rm = std::experimental::get< 3 >( energyRange );

    njoy::Log::info("\n Uranium-238"
                    "\n --------------"
                    "\n LRU: {}"
                    "\n LRF: {}"
                    "\n NRO: {}"
                    "\n NAPS: {}\n",
                    rm.LRU(), rm.LRF(), rm.NRO(), rm.NAPS() );

    auto& testData = std::get<1>( U238 );
    reichMoore::Apply{}( rm, test( testData ) );
  }
}

std::pair< njoy::ENDFtk::section::Type<2>, std::vector< double > >
resonances( const std::string& id ){
  auto testData = [&]{
    std::vector< double > data;
    std::ifstream tupleFile( id + "-tuple.txt" );

    if ( not tupleFile.is_open() ){
      njoy::Log::error( "Could not open file " + id + "-tuple.txt" );
      throw std::exception();
    }

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
