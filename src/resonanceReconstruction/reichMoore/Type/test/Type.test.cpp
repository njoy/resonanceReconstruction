#define CATCH_CONFIG_MAIN

#include <chrono>

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;
using namespace dimwits;

std::pair< njoy::ENDFtk::section::Type<2>, std::vector< double > >
resonances( const std::string& id );

static void escape( void* );
static void clobber();

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
    
    RANGES_FOR( auto energy, energies ){
      auto trial = xs( energy );
      escape( &trial );
    }
  }
}

auto timingTestExclusive( const std::vector< double >& testData ){
  return [ &testData ]( auto&& xs ){
    auto energies = testData | ranges::view::stride( 4 );
    
    auto start = std::chrono::high_resolution_clock::now();

    auto accumulation = 0.0*barns;
    RANGES_FOR( auto energy, energies ){
      accumulation += xs( energy );
    }
    escape( &accumulation );

    auto finish = std::chrono::high_resolution_clock::now();
    auto milliseconds =
      std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
    njoy::Log::info( "Approximately {} milliseconds passed while 'Apply'-ing",
                    milliseconds.count() );
  }
}

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
      Log::info( "inclusive timing" );
      auto start = std::chrono::high_resolution_clock::now();

      reichMoore::Apply{}( rm, timingTestInclusive( testData ) );

      auto finish = std::chrono::high_resolution_clock::now();
      auto milliseconds =
        std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
      njoy::Log::info( "Approximately {} milliseconds passed while 'Apply'-ing",
                      milliseconds.count() );
    }
    
    SECTION( "Exclusive" ){
      Log::info( "exclusive timing" );
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

    auto start = std::chrono::high_resolution_clock::now();

    reichMoore::Apply{}( rm, test( testData ) );

    auto finish = std::chrono::high_resolution_clock::now();
    auto milliseconds =
      std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
    njoy::Log::info( "Approximately {} milliseconds passed while 'Apply'-ing",
                    milliseconds.count() );
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

    auto start = std::chrono::high_resolution_clock::now();
    reichMoore::Apply{}( rm, test( testData ) );
    auto finish = std::chrono::high_resolution_clock::now();
    auto milliseconds =
      std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
    njoy::Log::info( "Approximately {} milliseconds passed while 'Apply'-ing",
                    milliseconds.count() );
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

    auto start = std::chrono::high_resolution_clock::now();
    reichMoore::Apply{}( rm, test( testData ) );
    auto finish = std::chrono::high_resolution_clock::now();
    auto milliseconds =
      std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
    njoy::Log::info( "Approximately {} milliseconds passed while 'Apply'-ing",
                    milliseconds.count() );
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

static void escape( void* p ){
  asm volatile( "" : : "g"(p) : "memory" );
}

static void clobber() {
  asm volatile( "" : : : "memory" );
}
