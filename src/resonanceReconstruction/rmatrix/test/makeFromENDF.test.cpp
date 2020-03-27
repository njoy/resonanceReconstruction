#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;
using namespace njoy::resonanceReconstruction::rmatrix;
using ParticlePairs = ENDF::resolved::RMatrixLimited::ParticlePairs;

constexpr AtomicMass neutronMass = 1.008664 * daltons;
constexpr ElectricalCharge elementaryCharge = 1.602e-19 * coulomb;

std::string chunkFromFe54();
std::string chunkFromCl35();

SCENARIO( "makeParticlePairs" ) {

  GIVEN( "a valid instance of ParticlePairs without errors" ) {

    std::string string = chunkFromFe54();
    auto begin = string.begin();
    auto end = string.end();
    long lineNumber = 1;

    ParticlePairs chunk( begin, end, lineNumber, 2625, 2, 151 );

    THEN( "a vector of ParticlePair is returned" ) {

      auto pairs = makeParticlePairs( chunk, neutronMass, elementaryCharge );

      CHECK( 2 == pairs.size() );

      const auto pair1 = pairs[0];
      CHECK( 0.0 == Approx( pair1.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair1.particle().charge().value ) );
      CHECK( 1.0 == Approx( pair1.particle().spin() ) );
      CHECK( +1 == pair1.particle().parity() );
      CHECK( 54.46635 * 1.008664 == Approx( pair1.residual().mass().value ) );
      CHECK( 26.0 * 1.602e-19 == Approx( pair1.residual().charge().value ) );
      CHECK( 0.0 == Approx( pair1.residual().spin() ) );
      CHECK( +1 == pair1.residual().parity() );
      CHECK( 0.0 == Approx( pair1.Q().value ) );
      CHECK( "mt102-a,mt102-b" == pair1.pairID() );

      const auto pair2 = pairs[1];
      CHECK( 1.008664 == Approx( pair2.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair2.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair2.particle().spin() ) );
      CHECK( +1 == pair2.particle().parity() );
      CHECK( 53.47624 * 1.008664 == Approx( pair2.residual().mass().value ) );
      CHECK( 26.0 * 1.602e-19 == Approx( pair2.residual().charge().value ) );
      CHECK( 0.0 == Approx( pair2.residual().spin() ) );
      CHECK( +1 == pair2.residual().parity() );
      CHECK( 0.0 == Approx( pair2.Q().value ) );
      CHECK( "mt2-a,mt2-b" == pair2.pairID() );
    } // THEN
  } // GIVEN
} // SCENARIO

std::string chunkFromFe54() {

  // particle pairs from Fe54 ENDF/B-VIII.0 LRF=7 resonance evaluation

  return
    " 0.000000+0 0.000000+0          2          0         24          42625 2151     \n"
    " 0.000000+0 5.446635+1 0.000000+0 2.600000+1 1.000000+0 0.000000+02625 2151     \n"
    " 0.000000+0 0.000000+0 0.000000+0 1.020000+2 0.000000+0 0.000000+02625 2151     \n"
    " 1.000000+0 5.347624+1 0.000000+0 2.600000+1 5.000000-1 0.000000+02625 2151     \n"
    " 0.000000+0 1.000000+0 0.000000+0 2.000000+0 0.000000+0 1.000000+02625 2151     \n";
}

std::string chunkFromCl35() {

  // particle pairs from Cl35 ENDF/B-VIII.0 LRF=7 resonance evaluation
  // particular features: Z for n,gamma is set to 0.0 (should be 17.0)

  return
    " 0.000000+0 0.000000+0          3          0         36          61725 2151     \n"
    " 0.000000+0 3.565932+1 0.000000+0 0.000000+0 1.000000+0 0.000000+01725 2151     \n"
    " 0.000000+0 0.000000+0 0.000000+0 1.020000+2 0.000000+0 0.000000+01725 2151     \n"
    " 1.000000+0 3.466845+1 0.000000+0 1.700000+1 5.000000-1 1.500000+01725 2151     \n"
    " 0.000000+0 1.000000+0 0.000000+0 2.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 9.986235-1 3.466863+1 1.000000+0 1.600000+1 5.000000-1 1.500000+01725 2151     \n"
    " 6.152200+5 1.000000+0 0.000000+0 6.000000+2 0.000000+0 0.000000+01725 2151     \n";
}
