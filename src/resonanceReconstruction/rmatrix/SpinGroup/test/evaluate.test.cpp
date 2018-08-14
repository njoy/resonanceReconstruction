#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using Particle = rmatrix::Particle;
using ParticlePair = rmatrix::ParticlePair;
using Neutron = rmatrix::Neutron;
using Photon = rmatrix::Photon;
using Fission = rmatrix::Fission;
template < typename Type > using Channel = rmatrix::Channel< Type >;
using Resonance = rmatrix::Resonance;
using ResonanceTable = rmatrix::ResonanceTable;
template < typename Option > using SpinGroup = rmatrix::SpinGroup< Option >;
using ReactionID = rmatrix::ReactionID;
using Sammy = rmatrix::Sammy;
using Constant = rmatrix::Constant;

constexpr AtomicMass neutronMass = 1.008664 * daltons;

SCENARIO( "evaluate" ) {

  //! @todo add test with more than one entrance channel
  //! @todo add test with a charged particle channel

  GIVEN( "valid data for a SpinGroup with one eliminated capture channel "
         "and one elastic channel" ) {

    // test based on Fe54 ENDF/B-VIII.0 LRF7 resonance evaluation
    // data given in Gamma = 2 gamma^2 P(Er) so conversion is required
    // cross section values extracted from NJOY2016.39

    // because the oribital angular momentum l = 0 for these SpinGroup,
    // SpinGroup< Sammy > and SpinGroup< Constant > should give the same results
    // as B = S = 0

    // using SpinGroup< Sammy > is equivalent to NJOY2016's LRF7 reconstruction

    // particles
    Particle photon( "g", 0.0 * daltons, 0.0 * coulombs, 1., +1);
    Particle neutron( "n", neutronMass, 0.0 * coulombs, 0.5, +1);
    Particle fe55( "Fe55_e0", 5.446635e+1 * neutronMass,
                              26.0 * coulombs, 0.0, +1);
    Particle fe54( "Fe54_e0", 5.347624e+1 * neutronMass,
                              26.0 * coulombs, 0.0, +1);

    // particle pairs
    ParticlePair in( neutron, fe54, 0.0 * electronVolt );
    ParticlePair out( photon, fe55, 0.0 * electronVolt );

    // channels
    Channel< Photon > capture( "1", out, { 0, 0.0, 0.5, +1 },
                               { 0.0 * rootBarn },
                               0.0 );
    Channel< Neutron > elastic( "2", in, { 0, 0.5, 0.5, +1 },
                                { 5.437300e-1 * rootBarn,
                                  5.437300e-1 * rootBarn },
                                0.0 );

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
    ResonanceTable single2 = single;

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
    ResonanceTable multiple2 = multiple;

    SpinGroup< Sammy > group1( in, { elastic }, std::move( single ) );
    SpinGroup< Sammy > group2( in, { elastic }, std::move( multiple ) );
    SpinGroup< Constant > group3( in, { elastic }, std::move( single2 ) );
    SpinGroup< Constant > group4( in, { elastic }, std::move( multiple2 ) );

    ReactionID elas = "n,Fe54_e0->n,Fe54_e0";
    ReactionID capt = "n,Fe54_e0->capture";

    THEN( "cross sections can be calculated for a single resonance using the "
          "Sammy boundary condition" ) {

      // first value is elastic, second value is eliminated capture
      tsl::hopscotch_map< ReactionID, Quantity< Barn > > xs;
      group1.evaluate( 1e-5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575880e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 6.895037e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e-4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575879e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.180402e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e-3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575878e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 6.895039e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e-2 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575861e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.180408e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e-1 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575695e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 6.895213e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+0 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.574032e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.180960e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+1 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.557416e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 6.912723e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+2 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.392161e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.237318e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.932783e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 9.067250e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 4.342056e+1 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.473549e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.979435e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 4.915666e-6 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+6 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.308817e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 1.343258e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 7.788000e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.424287e+2 == Approx( xs[ elas ].value ) );
      REQUIRE( 4.241421e-1 == Approx( xs[ capt ].value ) );
      xs.clear();
    }

    THEN( "cross sections can be calculated for multiple resonances using the "
          "Sammy boundary condition" ) {

      // first value is elastic, second value is eliminated capture
      tsl::hopscotch_map< ReactionID, Quantity< Barn > > xs;
      group2.evaluate( 1e-5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781791e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.082909e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e-4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781790e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.239813e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e-3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781781e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.082911e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e-2 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781682e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.239818e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e-1 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.780692e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.083086e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+0 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.770804e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.240372e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+1 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.672107e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.100642e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+2 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 7.704194e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.296876e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 7.098264e-3 == Approx( xs[ elas ].value ) );
      REQUIRE( 9.259109e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 4.062823e+1 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.502603e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 5.330137e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 5.716074e-5 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+6 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.324263e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 3.695041e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 7.788000e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.424287e+2 == Approx( xs[ elas ].value ) );
      REQUIRE( 4.241421e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 5.287200e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 4.738763e+1 == Approx( xs[ elas ].value ) );
      REQUIRE( 5.099758e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 7.190500e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.390972e+1 == Approx( xs[ elas ].value ) );
      REQUIRE( 4.208797e-2 == Approx( xs[ capt ].value ) );
      xs.clear();
    }

    THEN( "cross sections can be calculated for a single resonance using the "
          "constant boundary condition" ) {

      // first value is elastic, second value is eliminated capture
      tsl::hopscotch_map< ReactionID, Quantity< Barn > > xs;
      group3.evaluate( 1e-5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575880e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 6.895037e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e-4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575879e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.180402e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e-3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575878e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 6.895039e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e-2 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575861e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.180408e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e-1 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575695e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 6.895213e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+0 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.574032e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.180960e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+1 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.557416e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 6.912723e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+2 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.392161e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.237318e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.932783e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 9.067250e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 4.342056e+1 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.473549e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.979435e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 4.915666e-6 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+6 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.308817e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 1.343258e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 7.788000e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.424287e+2 == Approx( xs[ elas ].value ) );
      REQUIRE( 4.241421e-1 == Approx( xs[ capt ].value ) );
      xs.clear();
    }

    THEN( "cross sections can be calculated for multiple resonances using the "
          "constant boundary condition" ) {

      // first value is elastic, second value is eliminated capture
      tsl::hopscotch_map< ReactionID, Quantity< Barn > > xs;
      group4.evaluate( 1e-5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781791e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.082909e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e-4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781790e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.239813e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e-3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781781e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.082911e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e-2 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781682e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.239818e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e-1 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.780692e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.083086e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+0 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.770804e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.240372e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+1 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.672107e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.100642e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+2 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 7.704194e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.296876e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 7.098264e-3 == Approx( xs[ elas ].value ) );
      REQUIRE( 9.259109e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 4.062823e+1 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.502603e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 5.330137e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 5.716074e-5 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+6 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.324263e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 3.695041e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 7.788000e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.424287e+2 == Approx( xs[ elas ].value ) );
      REQUIRE( 4.241421e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 5.287200e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 4.738763e+1 == Approx( xs[ elas ].value ) );
      REQUIRE( 5.099758e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 7.190500e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.390972e+1 == Approx( xs[ elas ].value ) );
      REQUIRE( 4.208797e-2 == Approx( xs[ capt ].value ) );
      xs.clear();
    }
  } // GIVEN

  GIVEN( "valid data for a SpinGroup with one eliminated capture channel, "
         "one elastic channel and two fission channels" ) {

    // test based on Pu239 ENDF/B-VIII.0 LRF3 resonance evaluation (some of the
    // widths used in this test were negative in the original evaluation)
    // data given in Gamma = 2 gamma^2 P(Er) so conversion is required
    // cross section values extracted from NJOY2016.39 (note: LRF7 in NJOY2016
    // doesn't add potential scattering for missing J values)

    // because the oribital angular momentum l = 0 for these SpinGroup,
    // SpinGroup< Sammy > and SpinGroup< Constant > should give the same results

    // using SpinGroup< Sammy > is equivalent to NJOY2016

    // particles
    Particle photon( "g", 0.0 * daltons, 0.0 * coulombs, 1., +1);
    Particle neutron( "n", neutronMass, 0.0 * coulombs, 0.5, +1);
    Particle pu240( "Pu240_e0", 2.379916e+2 * neutronMass,
                                94.0 * coulombs, 0.5, +1);
    Particle pu239( "Pu239_e0", 2.369986e+2 * neutronMass,
                                94.0 * coulombs, 0.5, +1);
    Particle fission( "fission", 0.0 * daltons, 0.0 * coulombs, 0.0, +1);

    // particle pairs
    ParticlePair in( neutron, pu239, 0.0 * electronVolt );
    ParticlePair out1( photon, pu240, 0.0 * electronVolt );
    ParticlePair out2( fission, fission, 0.0 * electronVolt, "fission" );

    // channels
    Channel< Photon > capture( "1", out1, { 0, 0.0, 0.0, +1 },
                               { 0.0 * rootBarn },
                               0.0 );
    Channel< Neutron > elastic( "2", in, { 0, 0.5, 0.0, +1 },
                                { 9.410000e-1 * rootBarn },
                                0.0 );
    Channel< Fission > fission1( "3", out2, { 0, 0.0, 0.0, +1 },
                                 { 9.410000e-1 * rootBarn },
                                 0.0 );
    Channel< Fission > fission2( "4", out2, { 0, 0.0, 0.0, +1 },
                                 { 9.410000e-1 * rootBarn },
                                 0.0 );

    // conversion from Gamma to gamma
    auto eGamma = [&] ( double width, const Energy& energy ) -> ReducedWidth {
      return std::sqrt( width / 2. / elastic.penetrability( energy ) ) *
             rootElectronVolt;
    };
    auto fGamma = [&] ( double width ) -> ReducedWidth {
      return std::sqrt( width / 2. ) * rootElectronVolt;
    };
    auto cGamma = [&] ( double width ) -> ReducedWidth {
      return std::sqrt( width / 2. ) * rootElectronVolt;
    };

    // single resonance table
    ResonanceTable single(
      { "2", "3", "4" },
      { Resonance( 1.541700e+1 * electronVolt,
                   { eGamma( 2.056203e-3, 1.541700e+1 * electronVolt ),
                     fGamma( 1.093928e-6 ),
                     fGamma( 7.550000e-1 ) },
                   cGamma( 4.054259e-2 ) ) } );
    ResonanceTable single2 = single;

    // multiple resonance table
    ResonanceTable multiple(
      { "2", "3", "4" },
      { Resonance( 1.541700e+1 * electronVolt,
                   { eGamma( 2.056203e-3, 1.541700e+1 * electronVolt ),
                     fGamma( 1.093928e-6 ),
                     fGamma( 7.550000e-1 ) },
                   cGamma( 4.054259e-2 ) ),
        Resonance( 3.232700e+1 * electronVolt,
                   { eGamma( 8.678823e-4, 3.232700e+1 * electronVolt ),
                     fGamma( 5.235058e-3 ),
                     fGamma( 1.279000e-1 ) },
                   cGamma( 4.182541e-2 ) ),
        Resonance( 4.753400e+1 * electronVolt,
                   { eGamma( 5.171861e-3, 4.753400e+1 * electronVolt ),
                     fGamma( 5.548812e-1 ),
                     fGamma( 1.274000e-7 ) },
                   cGamma( 2.938826e-2 ) ) } );
    ResonanceTable multiple2 = multiple;

    SpinGroup< Sammy > group1( in, { elastic, fission1, fission2 }, 
                               std::move( single ) );
    SpinGroup< Sammy > group2( in, { elastic, fission1, fission2 }, 
                               std::move( multiple ) );
    SpinGroup< Constant > group3( in, { elastic, fission1, fission2 }, 
                                  std::move( single2 ) );
    SpinGroup< Constant > group4( in, { elastic, fission1, fission2 }, 
                                  std::move( multiple2 ) );

    ReactionID elas = "n,Pu239_e0->n,Pu239_e0";
    ReactionID fiss = "n,Pu239_e0->fission";
    ReactionID capt = "n,Pu239_e0->capture";

    THEN( "cross sections can be calculated for a single resonance using the "
          "Sammy boundary condition" ) {

      // first value is elastic, second and third value are fission and the
      // fourth value is eliminated capture
      tsl::hopscotch_map< ReactionID, Quantity< Barn > > xs;
      group1.evaluate( 1e-5 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.736137e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 8.627569e+1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 4.632894e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e-4 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.736137e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.728309e+1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 1.465067e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e-3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.736134e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 8.628676e+0 == Approx( xs[ fiss ].value ) );
      REQUIRE( 4.633488e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e-2 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.736107e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.731813e+0 == Approx( xs[ fiss ].value ) );
      REQUIRE( 1.466949e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e-1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.735839e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 8.740500e-1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 4.693537e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+0 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.732981e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 3.119579e-1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 1.675174e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.653310e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 6.955267e-1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 3.734889e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+2 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.789797e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 9.069863e-4 == Approx( xs[ fiss ].value ) );
      REQUIRE( 4.870400e-5 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.778610e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.116761e-6 == Approx( xs[ fiss ].value ) );
      REQUIRE( 1.136674e-7 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 2e+3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.774325e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 3.684029e-7 == Approx( xs[ fiss ].value ) );
      REQUIRE( 1.978276e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1.541700e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 3.050416e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 1.039122e+2 == Approx( xs[ fiss ].value ) );
      REQUIRE( 5.579952e+0 == Approx( xs[ capt ].value ) );
      xs.clear();
    }

    THEN( "cross sections can be calculated for multiple resonances using the "
          "Sammy boundary condition" ) {

      // first value is elastic, second and third value are fission and the
      // fourth value is eliminated capture
      tsl::hopscotch_map< ReactionID, Quantity< Barn > > xs;
      group2.evaluate( 1e-5 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.708738e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 1.164549e+2 == Approx( xs[ fiss ].value ) );
      REQUIRE( 5.455577e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e-4 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.708738e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 3.682667e+1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 1.725223e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e-3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.708734e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 1.164683e+1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 5.456212e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e-2 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.708702e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 3.686888e+0 == Approx( xs[ fiss ].value ) );
      REQUIRE( 1.727233e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e-1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.708378e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 1.178142e+0 == Approx( xs[ fiss ].value ) );
      REQUIRE( 5.520311e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+0 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.704945e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 4.151230e-1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 1.948730e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.618127e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.902586e-1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 3.880614e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+2 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.812236e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 4.139401e-3 == Approx( xs[ fiss ].value ) );
      REQUIRE( 2.030171e-4 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.779884e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 5.781515e-6 == Approx( xs[ fiss ].value ) );
      REQUIRE( 2.751723e-7 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 2e+3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.774946e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 9.897642e-7 == Approx( xs[ fiss ].value ) );
      REQUIRE( 4.705688e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1.541700e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 3.018352e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 1.039297e+2 == Approx( xs[ fiss ].value ) );
      REQUIRE( 5.580786e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 3.232700e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 3.200788e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.594132e+1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 2.384031e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 4.753400e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 3.826411e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 1.140900e+2 == Approx( xs[ fiss ].value ) );
      REQUIRE( 6.042426e+0 == Approx( xs[ capt ].value ) );
      xs.clear();
    }

    THEN( "cross sections can be calculated for a single resonance using the "
          "Sammy boundary condition" ) {

      // first value is elastic, second and third value are fission and the
      // fourth value is eliminated capture
      tsl::hopscotch_map< ReactionID, Quantity< Barn > > xs;
      group3.evaluate( 1e-5 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.736137e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 8.627569e+1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 4.632894e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e-4 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.736137e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.728309e+1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 1.465067e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e-3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.736134e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 8.628676e+0 == Approx( xs[ fiss ].value ) );
      REQUIRE( 4.633488e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e-2 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.736107e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.731813e+0 == Approx( xs[ fiss ].value ) );
      REQUIRE( 1.466949e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e-1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.735839e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 8.740500e-1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 4.693537e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+0 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.732981e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 3.119579e-1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 1.675174e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.653310e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 6.955267e-1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 3.734889e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+2 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.789797e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 9.069863e-4 == Approx( xs[ fiss ].value ) );
      REQUIRE( 4.870400e-5 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.778610e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.116761e-6 == Approx( xs[ fiss ].value ) );
      REQUIRE( 1.136674e-7 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 2e+3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.774325e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 3.684029e-7 == Approx( xs[ fiss ].value ) );
      REQUIRE( 1.978276e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1.541700e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 3.050416e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 1.039122e+2 == Approx( xs[ fiss ].value ) );
      REQUIRE( 5.579952e+0 == Approx( xs[ capt ].value ) );
      xs.clear();
    }

    THEN( "cross sections can be calculated for multiple resonances using the "
          "Sammy boundary condition" ) {

      // first value is elastic, second and third value are fission and the
      // fourth value is eliminated capture
      tsl::hopscotch_map< ReactionID, Quantity< Barn > > xs;
      group4.evaluate( 1e-5 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.708738e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 1.164549e+2 == Approx( xs[ fiss ].value ) );
      REQUIRE( 5.455577e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e-4 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.708738e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 3.682667e+1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 1.725223e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e-3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.708734e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 1.164683e+1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 5.456212e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e-2 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.708702e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 3.686888e+0 == Approx( xs[ fiss ].value ) );
      REQUIRE( 1.727233e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e-1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.708378e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 1.178142e+0 == Approx( xs[ fiss ].value ) );
      REQUIRE( 5.520311e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+0 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.704945e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 4.151230e-1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 1.948730e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.618127e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.902586e-1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 3.880614e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+2 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.812236e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 4.139401e-3 == Approx( xs[ fiss ].value ) );
      REQUIRE( 2.030171e-4 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.779884e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 5.781515e-6 == Approx( xs[ fiss ].value ) );
      REQUIRE( 2.751723e-7 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 2e+3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.774946e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 9.897642e-7 == Approx( xs[ fiss ].value ) );
      REQUIRE( 4.705688e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1.541700e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 3.018352e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 1.039297e+2 == Approx( xs[ fiss ].value ) );
      REQUIRE( 5.580786e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 3.232700e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 3.200788e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.594132e+1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 2.384031e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 4.753400e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 3.826411e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 1.140900e+2 == Approx( xs[ fiss ].value ) );
      REQUIRE( 6.042426e+0 == Approx( xs[ capt ].value ) );
      xs.clear();
    }
  } // GIVEN

  GIVEN( "valid data for a SpinGroup with a resonance at a negative energy" ) {

    // test based on Pu239 ENDF/B-VIII.0 LRF3 resonance evaluation (some of the
    // widths used in this test were negative in the original evaluation)
    // data given in Gamma = 2 gamma^2 P(Er) so conversion is required
    // cross section values extracted from NJOY2016.39 (note: LRF7 in NJOY2016
    // doesn't add potential scattering for missing J values)

    // particles
    Particle photon( "g", 0.0 * daltons, 0.0 * coulombs, 1., +1);
    Particle neutron( "n", neutronMass, 0.0 * coulombs, 0.5, +1);
    Particle pu240( "Pu240_e0", 2.379916e+2 * neutronMass,
                                94.0 * coulombs, 0.5, +1);
    Particle pu239( "Pu239_e0", 2.369986e+2 * neutronMass,
                                94.0 * coulombs, 0.5, +1);
    Particle fission( "fission", 0.0 * daltons, 0.0 * coulombs, 0.0, +1);

    // particle pairs
    ParticlePair in( neutron, pu239, 0.0 * electronVolt );
    ParticlePair out1( photon, pu240, 0.0 * electronVolt );
    ParticlePair out2( fission, fission, 0.0 * electronVolt, "fission" );

    // channels
    Channel< Photon > capture( "1", out1, { 0, 0.0, 0.0, +1 },
                               { 0.0 * rootBarn },
                               0.0 );
    Channel< Neutron > elastic( "2", in, { 0, 0.5, 0.0, +1 },
                                { 9.410000e-1 * rootBarn },
                                0.0 );
    Channel< Fission > fission1( "3", out2, { 0, 0.0, 0.0, +1 },
                                 { 9.410000e-1 * rootBarn },
                                 0.0 );
    Channel< Fission > fission2( "4", out2, { 0, 0.0, 0.0, +1 },
                                 { 9.410000e-1 * rootBarn },
                                 0.0 );

    // conversion from Gamma to gamma
    auto eGamma = [&] ( double width, const Energy& energy ) -> ReducedWidth {
      return std::sqrt( width / 2. / elastic.penetrability( energy ) ) *
             rootElectronVolt;
    };
    auto fGamma = [&] ( double width ) -> ReducedWidth {
      return std::sqrt( width / 2. ) * rootElectronVolt;
    };
    auto cGamma = [&] ( double width ) -> ReducedWidth {
      return std::sqrt( width / 2. ) * rootElectronVolt;
    };

    // multiple resonance table
    ResonanceTable table(
      { "2", "3", "4" },
      { Resonance( -1.541700e+1 * electronVolt,
                   { eGamma( 2.056203e-3, 1.541700e+1 * electronVolt ),
                     fGamma( 1.093928e-6 ),
                     fGamma( 7.550000e-1 ) },
                   cGamma( 4.054259e-2 ) ),
        Resonance( 3.232700e+1 * electronVolt,
                   { eGamma( 8.678823e-4, 3.232700e+1 * electronVolt ),
                     fGamma( 5.235058e-3 ),
                     fGamma( 1.279000e-1 ) },
                   cGamma( 4.182541e-2 ) ),
        Resonance( 4.753400e+1 * electronVolt,
                   { eGamma( 5.171861e-3, 4.753400e+1 * electronVolt ),
                     fGamma( 5.548812e-1 ),
                     fGamma( 1.274000e-7 ) },
                   cGamma( 2.938826e-2 ) ) } );

    SpinGroup< Sammy > group( in, { elastic, fission1, fission2 },
                              std::move( table ) );

    ReactionID elas = "n,Pu239_e0->n,Pu239_e0";
    ReactionID fiss = "n,Pu239_e0->fission";
    ReactionID capt = "n,Pu239_e0->capture";

    THEN( "cross sections can be calculated" ) {

      // first value is elastic, second and third value are fission and the
      // fourth value is eliminated capture
      tsl::hopscotch_map< ReactionID, Quantity< Barn > > xs;
      group.evaluate( 1e-5 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.800024e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.969822e+1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 5.456509e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e-4 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.800023e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.520250e+1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 1.725484e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e-3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.800020e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.968828e+0 == Approx( xs[ fiss ].value ) );
      REQUIRE( 5.455954e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e-2 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.799987e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.517112e+0 == Approx( xs[ fiss ].value ) );
      REQUIRE( 1.723733e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e-1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.799661e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.870530e-1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 5.401133e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e+0 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.796554e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.233478e-1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 1.566002e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.773379e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 3.532511e-2 == Approx( xs[ fiss ].value ) );
      REQUIRE( 3.181981e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e+2 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.809990e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 3.583571e-3 == Approx( xs[ fiss ].value ) );
      REQUIRE( 1.804712e-4 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e+3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.779863e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 5.625611e-6 == Approx( xs[ fiss ].value ) );
      REQUIRE( 2.683739e-7 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 2e+3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.774941e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 9.760312e-7 == Approx( xs[ fiss ].value ) );
      REQUIRE( 4.645619e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 3.232700e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 3.200741e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.593904e+1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 2.384228e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 4.753400e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 3.815638e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 1.140803e+2 == Approx( xs[ fiss ].value ) );
      REQUIRE( 6.042125e+0 == Approx( xs[ capt ].value ) );
      xs.clear();
    }
  } // GIVEN

  GIVEN( "valid data for a SpinGroup with resonances using negative widths" ) {

    // test based on Pu239 ENDF/B-VIII.0 LRF3 resonance evaluation (all
    // widths used in this test are from the original evaluation)
    // data given in Gamma = 2 gamma^2 P(Er) so conversion is required
    // cross section values extracted from NJOY2016.39 (note: LRF7 in NJOY2016
    // doesn't add potential scattering for missing J values)

    // particles
    Particle photon( "g", 0.0 * daltons, 0.0 * coulombs, 1., +1);
    Particle neutron( "n", neutronMass, 0.0 * coulombs, 0.5, +1);
    Particle pu240( "Pu240_e0", 2.379916e+2 * neutronMass,
                                94.0 * coulombs, 0.5, +1);
    Particle pu239( "Pu239_e0", 2.369986e+2 * neutronMass,
                                94.0 * coulombs, 0.5, +1);
    Particle fission( "fission", 0.0 * daltons, 0.0 * coulombs, 0.0, +1);

    // particle pairs
    ParticlePair in( neutron, pu239, 0.0 * electronVolt );
    ParticlePair out1( photon, pu240, 0.0 * electronVolt );
    ParticlePair out2( fission, fission, 0.0 * electronVolt, "fission" );

    // channels
    Channel< Photon > capture( "1", out1, { 0, 0.0, 0.0, +1 },
                               { 0.0 * rootBarn },
                               0.0 );
    Channel< Neutron > elastic( "2", in, { 0, 0.5, 0.0, +1 },
                                { 9.410000e-1 * rootBarn },
                                0.0 );
    Channel< Fission > fission1( "3", out2, { 0, 0.0, 0.0, +1 },
                                 { 9.410000e-1 * rootBarn },
                                 0.0 );
    Channel< Fission > fission2( "4", out2, { 0, 0.0, 0.0, +1 },
                                 { 9.410000e-1 * rootBarn },
                                 0.0 );

    // conversion from Gamma to gamma
    auto eGamma = [&] ( double width, const Energy& energy ) -> ReducedWidth {
      double sign = 1.0;
      if ( width < 0 ) sign = -1.0;
      return sign * std::sqrt( std::abs( width ) / 2. /
                               elastic.penetrability( energy ) ) *
             rootElectronVolt;
    };
    auto fGamma = [&] ( double width ) -> ReducedWidth {
      double sign = 1.0;
      if ( width < 0 ) sign = -1.0;
      return sign *std::sqrt( std::abs( width ) / 2. ) * rootElectronVolt;
    };
    auto cGamma = [&] ( double width ) -> ReducedWidth {
      return std::sqrt( width / 2. ) * rootElectronVolt;
    };

    // multiple resonance table
    ResonanceTable table(
      { "2", "3", "4" },
      { Resonance( 1.541700e+1 * electronVolt,
                   { eGamma( 2.056203e-3, 1.541700e+1 * electronVolt ),
                     fGamma( -1.093928e-6 ),
                     fGamma( 7.550000e-1 ) },
                   cGamma( 4.054259e-2 ) ),
        Resonance( 3.232700e+1 * electronVolt,
                   { eGamma( 8.678823e-4, 3.232700e+1 * electronVolt ),
                     fGamma( 5.235058e-3 ),
                     fGamma( -1.279000e-1 ) },
                   cGamma( 4.182541e-2 ) ),
        Resonance( 4.753400e+1 * electronVolt,
                   { eGamma( 5.171861e-3, 4.753400e+1 * electronVolt ),
                     fGamma( 5.548812e-1 ),
                     fGamma( -1.274000e-7 ) },
                   cGamma( 2.938826e-2 ) ) } );

    SpinGroup< Sammy > group( in, { elastic, fission1, fission2 },
                              std::move( table ) );

    ReactionID elas = "n,Pu239_e0->n,Pu239_e0";
    ReactionID fiss = "n,Pu239_e0->fission";
    ReactionID capt = "n,Pu239_e0->capture";

    THEN( "cross sections can be calculated" ) {

      // first value is elastic, second and third value are fission and the
      // fourth value is eliminated capture
      tsl::hopscotch_map< ReactionID, Quantity< Barn > > xs;
      group.evaluate( 1e-5 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.708724e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.968614e+1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 5.456366e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e-4 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.708723e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.519925e+1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 1.725473e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e-3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.708720e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.969599e+0 == Approx( xs[ fiss ].value ) );
      REQUIRE( 5.457001e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e-2 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.708688e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.523044e+0 == Approx( xs[ fiss ].value ) );
      REQUIRE( 1.727483e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e-1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.708364e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 8.069115e-1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 5.521115e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e+0 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.704930e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.868371e-1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 1.949035e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.617986e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 6.396225e-1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 3.881964e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e+2 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.812237e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 3.121466e-3 == Approx( xs[ fiss ].value ) );
      REQUIRE( 2.030192e-4 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e+3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.779884e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 3.851580e-6 == Approx( xs[ fiss ].value ) );
      REQUIRE( 2.751724e-7 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 2e+3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.774946e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 6.568293e-7 == Approx( xs[ fiss ].value ) );
      REQUIRE( 4.705688e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1.541700e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 2.983437e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 1.039338e+2 == Approx( xs[ fiss ].value ) );
      REQUIRE( 5.584260e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 3.232700e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 3.360875e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.591453e+1 == Approx( xs[ fiss ].value ) );
      REQUIRE( 2.389296e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 4.753400e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 3.826384e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 1.140768e+2 == Approx( xs[ fiss ].value ) );
      REQUIRE( 6.042619e+0 == Approx( xs[ capt ].value ) );
      xs.clear();
    }
  } // GIVEN
} // SCENARIO


