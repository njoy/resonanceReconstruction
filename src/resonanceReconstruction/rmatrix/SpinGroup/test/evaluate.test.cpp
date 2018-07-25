#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using Particle = rmatrix::Particle;
using ParticlePair = rmatrix::ParticlePair;
using Neutron = rmatrix::Neutron;
using Photon = rmatrix::Photon;
using Channel = rmatrix::Channel;
using Resonance = rmatrix::Resonance;
using ResonanceTable = rmatrix::ResonanceTable;
using SpinGroup = rmatrix::SpinGroup;

SCENARIO( "evaluate" ) {

  //! @todo add test with more than one channel (e.g. fission channels)
  //! @todo add test with a charged particle channel
  //! @todo add test with resonances at negative energies

  GIVEN( "valid data for a SpinGroup with one eliminated capture channel "
         "and one elastic channel" ) {

    // test based on Fe54 ENDF/B-VIII.0 LRF7 resonance evaluation
    // data given in Gamma = 2 gamma^2 P(Er) so conversion is required
    // cross section values extracted from NJOY2016.39

    // particles
    Particle photon( 0.0 * daltons, 0.0 * coulombs, 1., +1);
    Particle neutron( 1.008664 * daltons, 0.0 * coulombs, 0.5, +1);
    Particle fe55( 5.446635e+1 * 1.008664 * daltons, 26.0 * coulombs, 0.0, +1);
    Particle fe54( 5.347624e+1 * 1.008664 * daltons, 26.0 * coulombs, 0.0, +1);

    // particle pairs
    ParticlePair pair1( photon, fe55, 0.0 * electronVolt, "capture" );
    ParticlePair pair2( neutron, fe54, 0.0 * electronVolt, "elastic" );

    // channels
    Channel capture( "1", pair1, { 0, 0.0, 0.5, +1 },
                     { 0.0 * rootBarn },
                     0.0, Photon() );
    Channel elastic( "2", pair2, { 0, 0.5, 0.5, +1 },
                     { 5.437300e-1 * rootBarn, 5.437300e-1 * rootBarn },
                     0.0, Neutron() );

    // conversion from Gamma to gamma
    auto eGamma = [&] ( double width, const Energy& energy ) -> ReducedWidth {
      return std::sqrt( width / 2. / elastic.penetrability( energy ) ) *
             rootElectronVolt;
    };
    auto cGamma = [&] ( double width ) -> ReducedWidth {
      return std::sqrt( width / 2. ) * rootElectronVolt;
    };

    // single resonance table
    ResonanceTable single(
      { "2" },
      { Resonance( 7.788000e+3 * electronVolt,
                   { eGamma( 1.187354e+3, 7.788000e+3 * electronVolt ) },
                   cGamma( 1.455000e+0 ) ) } );

    // multiple resonance table
    ResonanceTable multiple(
      { "2" },
      { Resonance( 7.788000e+3 * electronVolt,
                   { eGamma( 1.187354e+3, 7.788000e+3 * electronVolt ) },
                   cGamma( 1.455000e+0 ) ),
        Resonance( 5.287200e+4 * electronVolt,
                   { eGamma( 2.000345e+3, 5.287200e+4 * electronVolt ) },
                   cGamma( 2.000000e+0 ) ),
        Resonance( 7.190500e+4 * electronVolt,
                   { eGamma( 1.781791e+3, 7.190500e+4 * electronVolt ) },
                   cGamma( 2.000000e+0 ) ) } );

    SpinGroup group1( { elastic }, std::move( single ) );
    SpinGroup group2( { elastic }, std::move( multiple ) );

    // simplifying gammas for tests
    auto getName = [] ( const auto& array, const unsigned int index )
                      { return std::get< 0 >( array[index] ); };
    auto getXS = [] ( const auto& array, const unsigned int index )
                    { return std::get< 1 >( array[index] ).value; };

    THEN( "cross sections can be calculated for a single resonance" ) {

      // first value is elastic, second value is eliminated capture
      auto xs = group1.evaluate( 1e-5 * electronVolt );
      REQUIRE( 2 == xs.size() );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 2.575880e-1 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 6.895037e+1 == Approx( getXS( xs, 1 ) ) );

      xs = group1.evaluate( 1e-4 * electronVolt );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 2.575879e-1 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 2.180402e+1 == Approx( getXS( xs, 1 ) ) );

      xs = group1.evaluate( 1e-3 * electronVolt );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 2.575878e-1 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 6.895039e+0 == Approx( getXS( xs, 1 ) ) );

      xs = group1.evaluate( 1e-2 * electronVolt );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 2.575861e-1 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 2.180408e+0 == Approx( getXS( xs, 1 ) ) );

      xs = group1.evaluate( 1e-1 * electronVolt );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 2.575695e-1 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 6.895213e-1 == Approx( getXS( xs, 1 ) ) );

      xs = group1.evaluate( 1e+0 * electronVolt );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 2.574032e-1 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 2.180960e-1 == Approx( getXS( xs, 1 ) ) );

      xs = group1.evaluate( 1e+1 * electronVolt );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 2.557416e-1 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 6.912723e-2 == Approx( getXS( xs, 1 ) ) );

      xs = group1.evaluate( 1e+2 * electronVolt );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 2.392161e-1 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 2.237318e-2 == Approx( getXS( xs, 1 ) ) );

      xs = group1.evaluate( 1e+3 * electronVolt );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 8.932783e-2 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 9.067250e-3 == Approx( getXS( xs, 1 ) ) );

      xs = group1.evaluate( 1e+4 * electronVolt );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 4.342056e+1 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 2.473549e-2 == Approx( getXS( xs, 1 ) ) );

      xs = group1.evaluate( 1e+5 * electronVolt );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 3.979435e+0 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 4.915666e-6 == Approx( getXS( xs, 1 ) ) );

      xs = group1.evaluate( 1e+6 * electronVolt );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 2.308817e+0 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 1.343258e-8 == Approx( getXS( xs, 1 ) ) );

      xs = group1.evaluate( 7.788000e+3 * electronVolt );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 3.424287e+2 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 4.241421e-1 == Approx( getXS( xs, 1 ) ) );
    }

    THEN( "cross sections can be calculated for multiple resonances" ) {

      // first value is elastic, second value is eliminated capture
      auto xs = group2.evaluate( 1e-5 * electronVolt );
      REQUIRE( 2 == xs.size() );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 8.781791e-2 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 7.082909e+1 == Approx( getXS( xs, 1 ) ) );

      xs = group2.evaluate( 1e-4 * electronVolt );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 8.781790e-2 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 2.239813e+1 == Approx( getXS( xs, 1 ) ) );

      xs = group2.evaluate( 1e-3 * electronVolt );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 8.781781e-2 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 7.082911e+0 == Approx( getXS( xs, 1 ) ) );

      xs = group2.evaluate( 1e-2 * electronVolt );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 8.781682e-2 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 2.239818e+0 == Approx( getXS( xs, 1 ) ) );

      xs = group2.evaluate( 1e-1 * electronVolt );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 8.780692e-2 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 7.083086e-1 == Approx( getXS( xs, 1 ) ) );

      xs = group2.evaluate( 1e+0 * electronVolt );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 8.770804e-2 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 2.240372e-1 == Approx( getXS( xs, 1 ) ) );

      xs = group2.evaluate( 1e+1 * electronVolt );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 8.672107e-2 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 7.100642e-2 == Approx( getXS( xs, 1 ) ) );

      xs = group2.evaluate( 1e+2 * electronVolt );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 7.704194e-2 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 2.296876e-2 == Approx( getXS( xs, 1 ) ) );

      xs = group2.evaluate( 1e+3 * electronVolt );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 7.098264e-3 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 9.259109e-3 == Approx( getXS( xs, 1 ) ) );

      xs = group2.evaluate( 1e+4 * electronVolt );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 4.062823e+1 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 2.502603e-2 == Approx( getXS( xs, 1 ) ) );

      xs = group2.evaluate( 1e+5 * electronVolt );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 5.330137e+0 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 5.716074e-5 == Approx( getXS( xs, 1 ) ) );

      xs = group2.evaluate( 1e+6 * electronVolt );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 2.324263e+0 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 3.695041e-8 == Approx( getXS( xs, 1 ) ) );

      xs = group2.evaluate( 7.788000e+3 * electronVolt );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 3.424287e+2 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 4.241421e-1 == Approx( getXS( xs, 1 ) ) );

      xs = group2.evaluate( 5.287200e+4 * electronVolt );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 4.738763e+1 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 5.099758e-2 == Approx( getXS( xs, 1 ) ) );

      xs = group2.evaluate( 7.190500e+4 * electronVolt );
      REQUIRE( "elastic" == getName( xs, 0 ) );
      REQUIRE( "capture" == getName( xs, 1 ) );
      REQUIRE( 3.390972e+1 == Approx( getXS( xs, 0 ) ) );
      REQUIRE( 4.208797e-2 == Approx( getXS( xs, 1 ) ) );
    }
  } // GIVEN
} // SCENARIO


