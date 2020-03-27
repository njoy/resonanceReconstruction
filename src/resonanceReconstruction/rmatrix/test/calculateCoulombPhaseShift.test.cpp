#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;
using namespace njoy::resonanceReconstruction::rmatrix;

// convenience typedefs
using OrbitalAngularMomentum = rmatrix::OrbitalAngularMomentum;

SCENARIO( "calculateCoulombPhaseShift" ) {

  GIVEN( "valid values for the orbital angular momentum and rho" ) {

    THEN( "the appropriate Coulomb phase shift is calculated" ) {

      // Neutron channels
      REQUIRE( 0.0 == Approx( calculateCoulombPhaseShift< Neutron >( 0, 0.5 ) ) );
      REQUIRE( 0.0 == Approx( calculateCoulombPhaseShift< Neutron >( 0, 1.0 ) ) );
      REQUIRE( 0.0 == Approx( calculateCoulombPhaseShift< Neutron >( 1, 0.5 ) ) );
      REQUIRE( 0.0 == Approx( calculateCoulombPhaseShift< Neutron >( 1, 1.0 ) ) );
      REQUIRE( 0.0 == Approx( calculateCoulombPhaseShift< Neutron >( 2, 0.5 ) ) );
      REQUIRE( 0.0 == Approx( calculateCoulombPhaseShift< Neutron >( 2, 1.0 ) ) );

      // Photon channels
      REQUIRE( 0.0 == Approx( calculateCoulombPhaseShift< Photon >( 0, 0.5 ) ) );
      REQUIRE( 0.0 == Approx( calculateCoulombPhaseShift< Photon >( 0, 1.0 ) ) );
      REQUIRE( 0.0 == Approx( calculateCoulombPhaseShift< Photon >( 1, 0.5 ) ) );
      REQUIRE( 0.0 == Approx( calculateCoulombPhaseShift< Photon >( 1, 1.0 ) ) );
      REQUIRE( 0.0 == Approx( calculateCoulombPhaseShift< Photon >( 2, 0.5 ) ) );
      REQUIRE( 0.0 == Approx( calculateCoulombPhaseShift< Photon >( 2, 1.0 ) ) );

      // Fission channels
      REQUIRE( 0.0 == Approx( calculateCoulombPhaseShift< Fission >( 0, 0.5 ) ) );
      REQUIRE( 0.0 == Approx( calculateCoulombPhaseShift< Fission >( 0, 1.0 ) ) );
      REQUIRE( 0.0 == Approx( calculateCoulombPhaseShift< Fission >( 1, 0.5 ) ) );
      REQUIRE( 0.0 == Approx( calculateCoulombPhaseShift< Fission >( 1, 1.0 ) ) );
      REQUIRE( 0.0 == Approx( calculateCoulombPhaseShift< Fission >( 2, 0.5 ) ) );
      REQUIRE( 0.0 == Approx( calculateCoulombPhaseShift< Fission >( 2, 1.0 ) ) );

      // ChargedParticle channels
      REQUIRE( 0.0 == Approx( calculateCoulombPhaseShift< ChargedParticle >( 0, 0.5 ) ) );
      REQUIRE( 0.0 == Approx( calculateCoulombPhaseShift< ChargedParticle >( 0, 1.0 ) ) );
      REQUIRE( 1.83048772171 == Approx( calculateCoulombPhaseShift< ChargedParticle >( 1, 0.5 ) ) );
      REQUIRE( 0.64209261593 == Approx( calculateCoulombPhaseShift< ChargedParticle >( 1, 1.0 ) ) );
      REQUIRE( 5.74680508636 == Approx( calculateCoulombPhaseShift< ChargedParticle >( 2, 0.5 ) ) );
      REQUIRE( 2.47258033765 == Approx( calculateCoulombPhaseShift< ChargedParticle >( 2, 1.0 ) ) );
    }
  } // GIVEN
} // SCENARIO
