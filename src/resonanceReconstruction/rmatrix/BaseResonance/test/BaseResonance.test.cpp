#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using BaseResonance = rmatrix::BaseResonance;

class TestResonance : public BaseResonance {

public:

  TestResonance( Energy energy, std::vector< ReducedWidth >&& widths ) :
    BaseResonance( energy, std::move( widths ) ) {}

  using BaseResonance::BaseResonance;
  using BaseResonance::energy;
  using BaseResonance::widths;
};

SCENARIO( "BaseResonance" ) {

  GIVEN( "valid data for a BaseResonance" ) {

    // single resonance data
    Energy energy = 6.823616e+4 * electronVolt;
    std::vector< ReducedWidth > widths = { 2.179040e+2 * rootElectronVolt,
                                           1.000000e-5 * rootElectronVolt };

    THEN( "a BaseResonance can be constructed" ) {

      TestResonance resonance( energy, std::move( widths ) );

      CHECK( 6.823616e+4 == Approx( resonance.energy().value ) );

      CHECK( 2 == resonance.widths().size() );
      CHECK( 2.179040e+2 == Approx( resonance.widths()[0].value ) );
      CHECK( 1.000000e-5 == Approx( resonance.widths()[1].value ) );
    } // THEN
  } // GIVEN
} // SCENARIO
