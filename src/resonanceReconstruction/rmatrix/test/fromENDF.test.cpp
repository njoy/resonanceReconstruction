#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;
using namespace njoy::resonanceReconstruction::rmatrix;
using ParticlePairs = ENDF::resolved::RMatrixLimited::ParticlePairs;
using ResonanceChannels = ENDF::resolved::RMatrixLimited::ResonanceChannels;

constexpr AtomicMass neutronMass = 1.008664 * daltons;
constexpr ElectricalCharge elementaryCharge = 1.602e-19 * coulomb;

std::string particlePairsFromFe54();
std::string resonanceChannelsFromFe54();
std::string particlePairsFromCl35();
std::string resonanceChannelsFromCl35();

SCENARIO( "makeParticlePairs" ) {

  GIVEN( "valid ENDF data" ) {

    std::string string = particlePairsFromFe54();
    auto begin = string.begin();
    auto end = string.end();
    long lineNumber = 1;

    ParticlePairs endfPairs( begin, end, lineNumber, 2625, 2, 151 );

    THEN( "a vector of ParticlePair is returned" ) {

      auto pairs = makeParticlePairs( endfPairs, neutronMass, elementaryCharge );

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

SCENARIO( "makeParticleChannels" ) {

  GIVEN( "valid ENDF data" ) {

    std::string string = particlePairsFromFe54();
    auto begin = string.begin();
    auto end = string.end();
    long lineNumber = 1;

    ParticlePairs endfPairs( begin, end, lineNumber, 2625, 2, 151 );
    std::vector< ParticlePair > pairs = makeParticlePairs( endfPairs, neutronMass, elementaryCharge );

    string = resonanceChannelsFromFe54();
    begin = string.begin();
    end = string.end();
    lineNumber = 1;

    ResonanceChannels endfChannels( begin, end, lineNumber, 2625, 2, 151 );

    THEN( "a vector of ParticleChannel is returned" ) {

      auto channels = makeParticleChannels( pairs, endfPairs, endfChannels );

      CHECK( 2 == channels.size() );

      const auto channel1 = std::get< Channel< Photon > >( channels[0] );
      const auto channel2 = std::get< Channel< Neutron > >( channels[1] );

      // particle pairs
      const auto pair1 = channel1.particlePair();
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

      const auto pair2 = channel2.particlePair();
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

      // quantum numbers
      const auto numbers1 = channel1.quantumNumbers();
      CHECK( 0 == numbers1.orbitalAngularMomentum() );
      CHECK( 0.0 == numbers1.spin() );
      CHECK( 0.5 == numbers1.totalAngularMomentum() );
      CHECK( +1 == numbers1.parity() );
      CHECK( "{0,0,1/2+}" == numbers1.toString() );

      const auto numbers2 = channel2.quantumNumbers();
      CHECK( 0 == numbers2.orbitalAngularMomentum() );
      CHECK( 0.5 == numbers2.spin() );
      CHECK( 0.5 == numbers2.totalAngularMomentum() );
      CHECK( +1 == numbers2.parity() );
      CHECK( "{0,1/2,1/2+}" == numbers2.toString() );

      // radii
      const auto radii1 = channel1.radii();
      CHECK( 0. == radii1.penetrabilityRadius( 1e-5 * electronVolt ).value );
      CHECK( 0. == radii1.shiftFactorRadius( 1e-5 * electronVolt ).value );
      CHECK( 0. == radii1.phaseShiftRadius( 1e-5 * electronVolt ).value );

      const auto radii2 = channel2.radii();
      CHECK( .54373 == Approx( radii2.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .54373 == Approx( radii2.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .54373 == Approx( radii2.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel1.boundaryCondition() );
      CHECK( 0. == channel2.boundaryCondition() );
    } // THEN
  } // GIVEN
} // SCENARIO

std::string particlePairsFromFe54() {

  // particle pairs from Fe54 ENDF/B-VIII.0 LRF=7 resonance evaluation

  return
    " 0.000000+0 0.000000+0          2          0         24          42625 2151     \n"
    " 0.000000+0 5.446635+1 0.000000+0 2.600000+1 1.000000+0 0.000000+02625 2151     \n"
    " 0.000000+0 0.000000+0 0.000000+0 1.020000+2 0.000000+0 0.000000+02625 2151     \n"
    " 1.000000+0 5.347624+1 0.000000+0 2.600000+1 5.000000-1 0.000000+02625 2151     \n"
    " 0.000000+0 1.000000+0 0.000000+0 2.000000+0 0.000000+0 1.000000+02625 2151     \n";
}

std::string resonanceChannelsFromFe54() {

  // resonance channels from Fe54 ENDF/B-VIII.0 LRF=7 resonance evaluation
  // J = 0.5 pi = +1

  return
    " 5.000000-1 0.000000+0          0          0         12          22625 2151     \n"
    " 1.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.000000+0 0.000000+0 5.000000-1 0.000000+0 5.437300-1 5.437300-12625 2151     \n";
}

std::string particlePairsFromCl35() {

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

std::string resonanceChannelsFromCl35() {

  // resonance channels from Cl35 ENDF/B-VIII.0 LRF=7 resonance evaluation
  // J = 1.0 pi = +1

  return

    " 1.000000+0 0.000000+0          0          0         18          31725 2151     \n"
    " 1.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.000000+0 0.000000+0 1.000000+0 0.000000+0 3.667980-1 4.822220-11725 2151     \n"
    " 3.000000+0 0.000000+0 1.000000+0 0.000000+0 3.667980-1 4.822220-11725 2151     \n";
}
