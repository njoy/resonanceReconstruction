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
  //! @todo add test with resonances at negative energies
  //! @todo add test with resonances with a negative width

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
    Particle photon( 0.0 * daltons, 0.0 * coulombs, 1., +1);
    Particle neutron( neutronMass, 0.0 * coulombs, 0.5, +1);
    Particle fe55( 5.446635e+1 * neutronMass, 26.0 * coulombs, 0.0, +1);
    Particle fe54( 5.347624e+1 * neutronMass, 26.0 * coulombs, 0.0, +1);

    // particle pairs
    ParticlePair pair1( photon, fe55, 0.0 * electronVolt, "capture" );
    ParticlePair pair2( neutron, fe54, 0.0 * electronVolt, "elastic" );

    // channels
    Channel< Photon > capture( "1", pair1, { 0, 0.0, 0.5, +1 },
                               { 0.0 * rootBarn },
                               0.0 );
    Channel< Neutron > elastic( "2", pair2, { 0, 0.5, 0.5, +1 },
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

    SpinGroup< Sammy > group1( { elastic }, { 0 }, std::move( single ) );
    SpinGroup< Sammy > group2( { elastic }, { 0 }, std::move( multiple ) );
    SpinGroup< Constant > group3( { elastic }, { 0 }, std::move( single2 ) );
    SpinGroup< Constant > group4( { elastic }, { 0 }, std::move( multiple2 ) );

    THEN( "cross sections can be calculated for a single resonance using the "
          "Sammy boundary condition" ) {

      // first value is elastic, second value is eliminated capture
      tsl::hopscotch_map< ReactionID, Quantity< Barn > > xs;
      group1.evaluate( 1e-5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575880e-1 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 6.895037e+1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group1.evaluate( 1e-4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575879e-1 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.180402e+1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group1.evaluate( 1e-3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575878e-1 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 6.895039e+0 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group1.evaluate( 1e-2 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575861e-1 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.180408e+0 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group1.evaluate( 1e-1 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575695e-1 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 6.895213e-1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group1.evaluate( 1e+0 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.574032e-1 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.180960e-1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group1.evaluate( 1e+1 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.557416e-1 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 6.912723e-2 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group1.evaluate( 1e+2 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.392161e-1 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.237318e-2 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group1.evaluate( 1e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.932783e-2 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 9.067250e-3 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group1.evaluate( 1e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 4.342056e+1 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.473549e-2 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group1.evaluate( 1e+5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.979435e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 4.915666e-6 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group1.evaluate( 1e+6 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.308817e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 1.343258e-8 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group1.evaluate( 7.788000e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.424287e+2 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 4.241421e-1 == Approx( xs[ "capture" ].value ) );
      xs.clear();
    }

    THEN( "cross sections can be calculated for multiple resonances using the "
          "Sammy boundary condition" ) {

      // first value is elastic, second value is eliminated capture
      tsl::hopscotch_map< ReactionID, Quantity< Barn > > xs;
      group2.evaluate( 1e-5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781791e-2 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 7.082909e+1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group2.evaluate( 1e-4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781790e-2 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.239813e+1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group2.evaluate( 1e-3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781781e-2 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 7.082911e+0 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group2.evaluate( 1e-2 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781682e-2 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.239818e+0 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group2.evaluate( 1e-1 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.780692e-2 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 7.083086e-1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group2.evaluate( 1e+0 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.770804e-2 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.240372e-1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group2.evaluate( 1e+1 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.672107e-2 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 7.100642e-2 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group2.evaluate( 1e+2 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 7.704194e-2 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.296876e-2 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group2.evaluate( 1e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 7.098264e-3 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 9.259109e-3 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group2.evaluate( 1e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 4.062823e+1 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.502603e-2 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group2.evaluate( 1e+5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 5.330137e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 5.716074e-5 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group2.evaluate( 1e+6 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.324263e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 3.695041e-8 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group2.evaluate( 7.788000e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.424287e+2 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 4.241421e-1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group2.evaluate( 5.287200e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 4.738763e+1 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 5.099758e-2 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group2.evaluate( 7.190500e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.390972e+1 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 4.208797e-2 == Approx( xs[ "capture" ].value ) );
      xs.clear();
    }

    THEN( "cross sections can be calculated for a single resonance using the "
          "constant boundary condition" ) {

      // first value is elastic, second value is eliminated capture
      tsl::hopscotch_map< ReactionID, Quantity< Barn > > xs;
      group3.evaluate( 1e-5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575880e-1 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 6.895037e+1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group3.evaluate( 1e-4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575879e-1 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.180402e+1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group3.evaluate( 1e-3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575878e-1 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 6.895039e+0 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group3.evaluate( 1e-2 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575861e-1 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.180408e+0 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group3.evaluate( 1e-1 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575695e-1 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 6.895213e-1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group3.evaluate( 1e+0 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.574032e-1 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.180960e-1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group3.evaluate( 1e+1 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.557416e-1 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 6.912723e-2 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group3.evaluate( 1e+2 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.392161e-1 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.237318e-2 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group3.evaluate( 1e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.932783e-2 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 9.067250e-3 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group3.evaluate( 1e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 4.342056e+1 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.473549e-2 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group3.evaluate( 1e+5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.979435e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 4.915666e-6 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group3.evaluate( 1e+6 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.308817e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 1.343258e-8 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group3.evaluate( 7.788000e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.424287e+2 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 4.241421e-1 == Approx( xs[ "capture" ].value ) );
      xs.clear();
    }

    THEN( "cross sections can be calculated for multiple resonances using the "
          "constant boundary condition" ) {

      // first value is elastic, second value is eliminated capture
      tsl::hopscotch_map< ReactionID, Quantity< Barn > > xs;
      group4.evaluate( 1e-5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781791e-2 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 7.082909e+1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group4.evaluate( 1e-4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781790e-2 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.239813e+1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group4.evaluate( 1e-3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781781e-2 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 7.082911e+0 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group4.evaluate( 1e-2 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781682e-2 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.239818e+0 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group4.evaluate( 1e-1 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.780692e-2 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 7.083086e-1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group4.evaluate( 1e+0 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.770804e-2 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.240372e-1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group4.evaluate( 1e+1 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.672107e-2 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 7.100642e-2 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group4.evaluate( 1e+2 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 7.704194e-2 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.296876e-2 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group4.evaluate( 1e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 7.098264e-3 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 9.259109e-3 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group4.evaluate( 1e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 4.062823e+1 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.502603e-2 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group4.evaluate( 1e+5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 5.330137e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 5.716074e-5 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group4.evaluate( 1e+6 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.324263e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 3.695041e-8 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group4.evaluate( 7.788000e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.424287e+2 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 4.241421e-1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group4.evaluate( 5.287200e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 4.738763e+1 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 5.099758e-2 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group4.evaluate( 7.190500e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.390972e+1 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 4.208797e-2 == Approx( xs[ "capture" ].value ) );
      xs.clear();
    }
  } // GIVEN

  GIVEN( "valid data for a SpinGroup with one eliminated capture channel, "
         "one elastic channel and two fission channels" ) {

    // test based on Pu239 ENDF/B-VIII.0 LRF3 resonance evaluation
    // data given in Gamma = 2 gamma^2 P(Er) so conversion is required
    // cross section values extracted from NJOY2016.39
    // the spin of Pu239 is set to 0.0 instead 0.5 to get a single J value test
    // (otherwise NJOY would add potential scattering for missing J values)

    // because the oribital angular momentum l = 0 for these SpinGroup,
    // SpinGroup< Sammy > and SpinGroup< Constant > should give the same results

    // using SpinGroup< Sammy > is equivalent to NJOY2016

    // particles
    Particle photon( 0.0 * daltons, 0.0 * coulombs, 1., +1);
    Particle neutron( neutronMass, 0.0 * coulombs, 0.5, +1);
    Particle pu240( 2.379916e+2 * neutronMass, 94.0 * coulombs, 0.0, +1);
    Particle pu239( 2.369986e+2 * neutronMass, 94.0 * coulombs, 0.0, +1);
    Particle fission( 0.0 * daltons, 0.0 * coulombs, 0.0, +1);

    // particle pairs
    ParticlePair pair1( photon, pu240, 0.0 * electronVolt, "capture" );
    ParticlePair pair2( neutron, pu239, 0.0 * electronVolt, "elastic" );
    ParticlePair pair3( fission, fission, 0.0 * electronVolt, "fission" );

    // channels
    Channel< Photon > capture( "1", pair1, { 0, 0.0, 0.0, +1 },
                               { 0.0 * rootBarn },
                               0.0 );
    Channel< Neutron > elastic( "2", pair2, { 0, 0.5, 0.0, +1 },
                                { 9.410000e-1 * rootBarn },
                                0.0 );
    Channel< Fission > fission1( "3", pair3, { 0, 0.0, 0.0, +1 },
                                 { 9.410000e-1 * rootBarn },
                                 0.0 );
    Channel< Fission > fission2( "4", pair3, { 0, 0.0, 0.0, +1 },
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

    SpinGroup< Sammy > group1( { elastic, fission1, fission2 }, { 0 },
                               std::move( single ) );
    SpinGroup< Sammy > group2( { elastic, fission1, fission2 }, { 0 },
                               std::move( multiple ) );
    SpinGroup< Constant > group3( { elastic, fission1, fission2 }, { 0 },
                                  std::move( single2 ) );
    SpinGroup< Constant > group4( { elastic, fission1, fission2 }, { 0 },
                                 std::move( multiple2 ) );

    THEN( "cross sections can be calculated for a single resonance using the "
          "Sammy boundary condition" ) {

      // first value is elastic, second and third value are fission and the
      // fourth value is eliminated capture
      tsl::hopscotch_map< ReactionID, Quantity< Barn > > xs;
      group1.evaluate( 1e-5 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.472274e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 1.725514e+2 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 9.265788e+0 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group1.evaluate( 1e-4 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.472274e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 5.456617e+1 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 2.930134e+0 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group1.evaluate( 1e-3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.472268e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 1.725735e+1 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 9.266977e-1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group1.evaluate( 1e-2 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.472215e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 5.463627e+0 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 2.933898e-1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group1.evaluate( 1e-1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.471679e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 1.748100e+0 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 9.387073e-2 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group1.evaluate( 1e+0 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.465961e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 6.239157e-1 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 3.350347e-2 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group1.evaluate( 1e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.306619e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 1.391053e+0 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 7.469779e-2 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group1.evaluate( 1e+2 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.579595e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 1.813973e-3 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 9.740800e-5 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group1.evaluate( 1e+3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.557221e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 4.233523e-6 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 2.273347e-7 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group1.evaluate( 2e+3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.548651e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 7.368059e-7 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 3.956553e-8 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group1.evaluate( 1.541700e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 6.100833e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.078244e+2 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 1.115990e+1 == Approx( xs[ "capture" ].value ) );
      xs.clear();
    }

    THEN( "cross sections can be calculated for multiple resonances using the "
          "Sammy boundary condition" ) {

      // first value is elastic, second and third value are fission and the
      // fourth value is eliminated capture
      tsl::hopscotch_map< ReactionID, Quantity< Barn > > xs;
      group2.evaluate( 1e-5 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.417476e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.329099e+2 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 1.091115e+1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group2.evaluate( 1e-4 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.417475e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 7.365334e+1 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 3.450446e+0 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group2.evaluate( 1e-3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.417469e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.329366e+1 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 1.091242e+0 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group2.evaluate( 1e-2 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.417404e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 7.373775e+0 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 3.454466e-1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group2.evaluate( 1e-1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.416755e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.356285e+0 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 1.104062e-1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group2.evaluate( 1e+0 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.409891e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 8.302460e-1 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 3.897461e-2 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group2.evaluate( 1e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.236255e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 1.580517e+0 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 7.761229e-2 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group2.evaluate( 1e+2 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.624471e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 8.278803e-3 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 4.060342e-4 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group2.evaluate( 1e+3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.559769e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 1.156303e-5 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 5.503447e-7 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group2.evaluate( 2e+3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.549892e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 1.979528e-6 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 9.411375e-8 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group2.evaluate( 1.541700e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 6.036703e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.078595e+2 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 1.116157e+1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group2.evaluate( 3.232700e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 6.401576e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 1.518826e+2 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 4.768063e+1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group2.evaluate( 4.753400e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 7.652821e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.281801e+2 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 1.208485e+1 == Approx( xs[ "capture" ].value ) );
      xs.clear();
    }

    THEN( "cross sections can be calculated for a single resonance using the "
          "Sammy boundary condition" ) {

      // first value is elastic, second and third value are fission and the
      // fourth value is eliminated capture
      tsl::hopscotch_map< ReactionID, Quantity< Barn > > xs;
      group3.evaluate( 1e-5 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.472274e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 1.725514e+2 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 9.265788e+0 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group3.evaluate( 1e-4 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.472274e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 5.456617e+1 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 2.930134e+0 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group3.evaluate( 1e-3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.472268e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 1.725735e+1 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 9.266977e-1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group3.evaluate( 1e-2 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.472215e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 5.463627e+0 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 2.933898e-1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group3.evaluate( 1e-1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.471679e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 1.748100e+0 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 9.387073e-2 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group3.evaluate( 1e+0 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.465961e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 6.239157e-1 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 3.350347e-2 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group3.evaluate( 1e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.306619e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 1.391053e+0 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 7.469779e-2 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group3.evaluate( 1e+2 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.579595e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 1.813973e-3 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 9.740800e-5 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group3.evaluate( 1e+3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.557221e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 4.233523e-6 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 2.273347e-7 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group3.evaluate( 2e+3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.548651e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 7.368059e-7 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 3.956553e-8 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group3.evaluate( 1.541700e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 6.100833e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.078244e+2 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 1.115990e+1 == Approx( xs[ "capture" ].value ) );
      xs.clear();
    }

    THEN( "cross sections can be calculated for multiple resonances using the "
          "Sammy boundary condition" ) {

      // first value is elastic, second and third value are fission and the
      // fourth value is eliminated capture
      tsl::hopscotch_map< ReactionID, Quantity< Barn > > xs;
      group4.evaluate( 1e-5 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.417476e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.329099e+2 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 1.091115e+1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group4.evaluate( 1e-4 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.417475e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 7.365334e+1 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 3.450446e+0 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group4.evaluate( 1e-3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.417469e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.329366e+1 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 1.091242e+0 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group4.evaluate( 1e-2 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.417404e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 7.373775e+0 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 3.454466e-1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group4.evaluate( 1e-1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.416755e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.356285e+0 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 1.104062e-1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group4.evaluate( 1e+0 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.409891e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 8.302460e-1 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 3.897461e-2 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group4.evaluate( 1e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.236255e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 1.580517e+0 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 7.761229e-2 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group4.evaluate( 1e+2 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.624471e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 8.278803e-3 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 4.060342e-4 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group4.evaluate( 1e+3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.559769e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 1.156303e-5 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 5.503447e-7 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group4.evaluate( 2e+3 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 5.549892e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 1.979528e-6 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 9.411375e-8 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group4.evaluate( 1.541700e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 6.036703e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.078595e+2 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 1.116157e+1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group4.evaluate( 3.232700e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 6.401576e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 1.518826e+2 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 4.768063e+1 == Approx( xs[ "capture" ].value ) );
      xs.clear();

      group4.evaluate( 4.753400e+1 * electronVolt, xs );
      REQUIRE( 3 == xs.size() );
      REQUIRE( 7.652821e+0 == Approx( xs[ "elastic" ].value ) );
      REQUIRE( 2.281801e+2 == Approx( xs[ "fission" ].value ) );
      REQUIRE( 1.208485e+1 == Approx( xs[ "capture" ].value ) );
      xs.clear();
    }
  } // GIVEN
} // SCENARIO


