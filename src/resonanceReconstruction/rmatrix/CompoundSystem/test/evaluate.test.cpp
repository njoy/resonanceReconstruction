#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using Particle = rmatrix::Particle;
using ParticlePair = rmatrix::ParticlePair;
using Neutron = rmatrix::Neutron;
using Photon = rmatrix::Photon;
using Fission = rmatrix::Fission;
using ReactionID = rmatrix::ReactionID;
template < typename Type > using Channel = rmatrix::Channel< Type >;
using Resonance = rmatrix::Resonance;
using ResonanceTable = rmatrix::ResonanceTable;
using Sammy = rmatrix::Sammy;
using Constant = rmatrix::Constant;
template < typename Option > using SpinGroup = rmatrix::SpinGroup< Option >;
template < typename Option > using CompoundSystem = rmatrix::CompoundSystem< Option >;

constexpr AtomicMass neutronMass = 1.008664 * daltons;

SCENARIO( "evaluate" ) {

  GIVEN( "valid data for a CompoundSystem with only one SpinGroup without "
         "missing J values" ) {

    // test based on Fe54 ENDF/B-VIII.0 LRF7 resonance evaluation
    // data given in Gamma = 2 gamma^2 P(Er) so conversion is required
    // cross section values extracted from NJOY2016.39

    // because the oribital angular momentum l = 0 for these SpinGroup,
    // CompoundSystem< Sammy > and CompoundSystem< Constant > should give the
    // same results as B = S = 0

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
    Channel< Photon > capture( out, { 0, 0.0, 0.5, +1 },
                               { 0.0 * rootBarn },
                               0.0 );
    Channel< Neutron > elastic( in, { 0, 0.5, 0.5, +1 },
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
      { elastic.channelID() },
      { Resonance( 7.788000e+3 * electronVolt,
                   { eGamma( 1.187354e+3, 7.788000e+3 * electronVolt ) },
                   cGamma( 1.455000e+0 ) ) } );
    ResonanceTable single2 = single;

    // multiple resonance table
    ResonanceTable multiple(
      { elastic.channelID() },
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

    CompoundSystem< Sammy > system1( { group1 } );
    CompoundSystem< Sammy > system2( { group2 } );
    CompoundSystem< Constant > system3( { group3 } );
    CompoundSystem< Constant > system4( { group4 } );

    ReactionID elas = "n,Fe54_e0->n,Fe54_e0";
    ReactionID capt = "n,Fe54_e0->capture";

    THEN( "cross sections can be calculated for a single resonance using the "
          "Sammy boundary condition" ) {

      // first value is elastic, second value is eliminated capture
      tsl::hopscotch_map< ReactionID, Quantity< Barn > > xs;
      system1.evaluate( 1e-5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575880e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 6.895037e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system1.evaluate( 1e-4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575879e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.180402e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system1.evaluate( 1e-3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575878e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 6.895039e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      system1.evaluate( 1e-2 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575861e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.180408e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      system1.evaluate( 1e-1 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575695e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 6.895213e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system1.evaluate( 1e+0 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.574032e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.180960e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system1.evaluate( 1e+1 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.557416e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 6.912723e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system1.evaluate( 1e+2 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.392161e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.237318e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system1.evaluate( 1e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.932783e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 9.067250e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      system1.evaluate( 1e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 4.342056e+1 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.473549e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system1.evaluate( 1e+5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.979435e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 4.915666e-6 == Approx( xs[ capt ].value ) );
      xs.clear();

      system1.evaluate( 1e+6 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.308817e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 1.343258e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      system1.evaluate( 7.788000e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.424287e+2 == Approx( xs[ elas ].value ) );
      REQUIRE( 4.241421e-1 == Approx( xs[ capt ].value ) );
      xs.clear();
    }

    THEN( "cross sections can be calculated for multiple resonances using the "
          "Sammy boundary condition" ) {

      // first value is elastic, second value is eliminated capture
      tsl::hopscotch_map< ReactionID, Quantity< Barn > > xs;
      system2.evaluate( 1e-5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781791e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.082909e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 1e-4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781790e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.239813e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 1e-3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781781e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.082911e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 1e-2 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781682e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.239818e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 1e-1 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.780692e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.083086e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 1e+0 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.770804e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.240372e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 1e+1 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.672107e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.100642e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 1e+2 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 7.704194e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.296876e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 1e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 7.098264e-3 == Approx( xs[ elas ].value ) );
      REQUIRE( 9.259109e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 1e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 4.062823e+1 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.502603e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 1e+5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 5.330137e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 5.716074e-5 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 1e+6 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.324263e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 3.695041e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 7.788000e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.424287e+2 == Approx( xs[ elas ].value ) );
      REQUIRE( 4.241421e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 5.287200e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 4.738763e+1 == Approx( xs[ elas ].value ) );
      REQUIRE( 5.099758e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 7.190500e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.390972e+1 == Approx( xs[ elas ].value ) );
      REQUIRE( 4.208797e-2 == Approx( xs[ capt ].value ) );
      xs.clear();
    }

    THEN( "cross sections can be calculated for a single resonance using the "
          "constant boundary condition" ) {

      // first value is elastic, second value is eliminated capture
      tsl::hopscotch_map< ReactionID, Quantity< Barn > > xs;
      system3.evaluate( 1e-5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575880e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 6.895037e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system3.evaluate( 1e-4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575879e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.180402e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system3.evaluate( 1e-3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575878e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 6.895039e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      system3.evaluate( 1e-2 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575861e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.180408e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      system3.evaluate( 1e-1 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.575695e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 6.895213e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system3.evaluate( 1e+0 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.574032e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.180960e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system3.evaluate( 1e+1 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.557416e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 6.912723e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system3.evaluate( 1e+2 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.392161e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.237318e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system3.evaluate( 1e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.932783e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 9.067250e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      system3.evaluate( 1e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 4.342056e+1 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.473549e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system3.evaluate( 1e+5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.979435e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 4.915666e-6 == Approx( xs[ capt ].value ) );
      xs.clear();

      system3.evaluate( 1e+6 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.308817e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 1.343258e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      system3.evaluate( 7.788000e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.424287e+2 == Approx( xs[ elas ].value ) );
      REQUIRE( 4.241421e-1 == Approx( xs[ capt ].value ) );
      xs.clear();
    }

    THEN( "cross sections can be calculated for multiple resonances using the "
          "constant boundary condition" ) {

      // first value is elastic, second value is eliminated capture
      tsl::hopscotch_map< ReactionID, Quantity< Barn > > xs;
      system4.evaluate( 1e-5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781791e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.082909e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 1e-4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781790e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.239813e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 1e-3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781781e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.082911e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 1e-2 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781682e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.239818e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 1e-1 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.780692e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.083086e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 1e+0 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.770804e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.240372e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 1e+1 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.672107e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.100642e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 1e+2 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 7.704194e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.296876e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 1e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 7.098264e-3 == Approx( xs[ elas ].value ) );
      REQUIRE( 9.259109e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 1e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 4.062823e+1 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.502603e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 1e+5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 5.330137e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 5.716074e-5 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 1e+6 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.324263e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 3.695041e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 7.788000e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.424287e+2 == Approx( xs[ elas ].value ) );
      REQUIRE( 4.241421e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 5.287200e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 4.738763e+1 == Approx( xs[ elas ].value ) );
      REQUIRE( 5.099758e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 7.190500e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.390972e+1 == Approx( xs[ elas ].value ) );
      REQUIRE( 4.208797e-2 == Approx( xs[ capt ].value ) );
      xs.clear();
    }
  } // GIVEN

  GIVEN( "valid data for a CompoundSystem with five SpinGroup without "
         "missing J values" ) {

    // test based on Fe54 ENDF/B-VIII.0 LRF7 resonance evaluation
    // data given in Gamma = 2 gamma^2 P(Er) so conversion is required
    // cross section values extracted from NJOY2016.39

    // because the oribital angular momentum values for these SpinGroup are
    // 0, 1 and 2, CompoundSystem< Sammy > and CompoundSystem< Constant > will
    // give different results

    // using SpinGroup< Sammy > is equivalent to NJOY2016's LRF7 reconstruction

    // particles
    Particle photon( "g", 0.0 * daltons, 0.0 * coulombs, 1., +1);
    Particle neutron( "n", neutronMass, 0.0 * coulombs, 0.5, +1);
    Particle fe55( "Fe55_e0", 5.446635e+1 * neutronMass,
                              26.0 * coulombs, 0.0, +1);
    Particle fe54( "Fe54_e0", 5.347624e+1 * neutronMass,
                              26.0 * coulombs, 0.0, +1);

    // particle pairs
    ParticlePair out( photon, fe55, 0.0 * electronVolt );
    ParticlePair in( neutron, fe54, 0.0 * electronVolt );

    // channels
    Channel< Photon > capture1( out, { 0, 0.0, 0.5, +1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic1( in, { 0, 0.5, 0.5, +1 },
                                 { 5.437300e-1 * rootBarn,
                                   5.437300e-1 * rootBarn } );
    Channel< Photon > capture2( out, { 0, 0.0, 0.5, -1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic2( in, { 1, 0.5, 0.5, -1 },
                                 { 5.437300e-1 * rootBarn,
                                   5.437300e-1 * rootBarn } );
    Channel< Photon > capture3( out, { 0, 0.0, 1.5, -1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic3( in, { 1, 0.5, 1.5, -1 },
                                 { 5.437300e-1 * rootBarn,
                                   5.437300e-1 * rootBarn } );
    Channel< Photon > capture4( out, { 0, 0.0, 1.5, +1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic4( in, { 2, 0.5, 1.5, +1 },
                                 { 5.437300e-1 * rootBarn,
                                   5.437300e-1 * rootBarn } );
    Channel< Photon > capture5( out, { 0, 0.0, 2.5, +1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic5( in, { 2, 0.5, 2.5, +1 },
                                 { 5.437300e-1 * rootBarn,
                                   5.437300e-1 * rootBarn } );

    // conversion from Gamma to gamma
    auto eGamma = [&] ( double width, const Energy& energy,
                        const Channel< Neutron >& elastic ) -> ReducedWidth {
      return std::sqrt( width / 2. / elastic.penetrability( energy ) ) *
             rootElectronVolt;
    };
    auto cGamma = [&] ( double width ) -> ReducedWidth {
      return std::sqrt( width / 2. ) * rootElectronVolt;
    };

    // multiple resonance table
    ResonanceTable table1(
      { elastic1.channelID() },
      { Resonance( 7.788000e+3 * electronVolt,
                   { eGamma( 1.187354e+3, 7.788000e+3 * electronVolt,
                             elastic1 ) },
                   cGamma( 1.455000e+0 ) ),
        Resonance( 5.287200e+4 * electronVolt,
                   { eGamma( 2.000345e+3, 5.287200e+4 * electronVolt,
                             elastic1 ) },
                   cGamma( 2.000000e+0 ) ),
        Resonance( 7.190500e+4 * electronVolt,
                   { eGamma( 1.781791e+3, 7.190500e+4 * electronVolt,
                             elastic1 ) },
                   cGamma( 2.000000e+0 ) ) } );
    ResonanceTable table2(
      { elastic2.channelID() },
      { Resonance( 5.152000e+4 * electronVolt,
                   { eGamma( 1.600200e+1, 5.152000e+4 * electronVolt,
                             elastic2 ) },
                   cGamma( 3.600000e-1 ) ),
        Resonance( 5.359000e+4 * electronVolt,
                   { eGamma( 1.700000e+1, 5.359000e+4 * electronVolt,
                             elastic2 ) },
                   cGamma( 1.500000e+0 ) ),
        Resonance( 5.545900e+4 * electronVolt,
                   { eGamma( 3.200000e+1, 5.545900e+4 * electronVolt,
                             elastic2 ) },
                   cGamma( 5.600000e-1 ) ) } );
    ResonanceTable table3(
      { elastic3.channelID() },
      { Resonance( 3.099000e+3 * electronVolt,
                   { eGamma( 1.400000e-3, 3.099000e+3 * electronVolt,
                             elastic3 ) },
                   cGamma( 5.900000e-1 ) ),
        Resonance( 1.358100e+4 * electronVolt,
                   { eGamma( 1.750000e-2, 1.358100e+4 * electronVolt,
                             elastic3 ) },
                   cGamma( 5.900000e-1 ) ),
        Resonance( 1.927800e+4 * electronVolt,
                   { eGamma( 2.750000e-2, 1.927800e+4 * electronVolt,
                             elastic3 ) },
                   cGamma( 5.900000e-1 ) ) } );
    ResonanceTable table4(
      { elastic4.channelID() },
      { Resonance( 9.480000e+3 * electronVolt,
                   { eGamma( 1.200000e+0, 9.480000e+3 * electronVolt,
                             elastic4 ) },
                   cGamma( 2.700000e-1 ) ),
        Resonance( 1.118000e+4 * electronVolt,
                   { eGamma( 3.850100e+0, 1.118000e+4 * electronVolt,
                             elastic4 ) },
                   cGamma( 3.500000e-1 ) ),
        Resonance( 1.445000e+4 * electronVolt,
                   { eGamma( 7.000200e-1, 1.445000e+4 * electronVolt,
                             elastic4 ) },
                   cGamma( 3.500000e-1 ) ) } );
    ResonanceTable table5(
      { elastic5.channelID() },
      { Resonance( 1.264000e+5 * electronVolt,
                   { eGamma( 2.900000e+0, 1.264000e+5 * electronVolt,
                             elastic5 ) },
                   cGamma( 1.100000e+0 ) ),
        Resonance( 1.504700e+5 * electronVolt,
                   { eGamma( 2.600000e+0, 1.504700e+5 * electronVolt,
                             elastic5 ) },
                   cGamma( 9.600000e-1 ) ),
        Resonance( 1.779400e+5 * electronVolt,
                   { eGamma( 1.400000e+0, 1.779400e+5 * electronVolt,
                             elastic5 ) },
                   cGamma( 9.600000e-1 ) ) } );

    SpinGroup< Sammy > group1( in, { elastic1 }, std::move( table1 ) );
    SpinGroup< Sammy > group2( in, { elastic2 }, std::move( table2 ) );
    SpinGroup< Sammy > group3( in, { elastic3 }, std::move( table3 ) );
    SpinGroup< Sammy > group4( in, { elastic4 }, std::move( table4 ) );
    SpinGroup< Sammy > group5( in, { elastic5 }, std::move( table5 ) );

    CompoundSystem< Sammy > system( { group1, group2, group3, group4, group5 } );

    ReactionID elas = "n,Fe54_e0->n,Fe54_e0";
    ReactionID capt = "n,Fe54_e0->capture";

    THEN( "cross sections can be calculated" ) {

      // first value is elastic, second value is eliminated capture
      tsl::hopscotch_map< ReactionID, Quantity< Barn > > xs;
      system.evaluate( 1e-5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781791e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.082909e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e-4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781790e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.239813e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e-3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781781e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.082911e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e-2 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.781682e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.239818e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e-1 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.780692e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.083086e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e+0 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.770804e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.240372e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e+1 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.672107e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 7.100642e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e+2 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 7.704196e-2 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.296878e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 7.100443e-3 == Approx( xs[ elas ].value ) );
      REQUIRE( 9.259257e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 4.062844e+1 == Approx( xs[ elas ].value ) );
      REQUIRE( 2.531089e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e+5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 5.355058e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 5.836614e-5 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e+6 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.240021e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 1.613911e-7 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 7.788000e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.424289e+2 == Approx( xs[ elas ].value ) );
      REQUIRE( 4.241627e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 5.287200e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 4.739413e+1 == Approx( xs[ elas ].value ) );
      REQUIRE( 5.169252e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 7.190500e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.392478e+1 == Approx( xs[ elas ].value ) );
      REQUIRE( 4.209119e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 5.152000e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 5.688152e+1 == Approx( xs[ elas ].value ) );
      REQUIRE( 1.147127e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 5.359000e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 8.699519e+1 == Approx( xs[ elas ].value ) );
      REQUIRE( 3.790368e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 5.545900e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 6.447107e+1 == Approx( xs[ elas ].value ) );
      REQUIRE( 8.302520e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 3.099000e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 4.288507e-1 == Approx( xs[ elas ].value ) );
      REQUIRE( 4.129389e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1.358100e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 1.269370e+1 == Approx( xs[ elas ].value ) );
      REQUIRE( 1.113654e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1.927800e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 7.072284e+0 == Approx( xs[ elas ].value ) );
      REQUIRE( 1.192963e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 9.480000e+3 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 4.384506e+2 == Approx( xs[ elas ].value ) );
      REQUIRE( 8.551893e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1.118000e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 4.293224e+2 == Approx( xs[ elas ].value ) );
      REQUIRE( 3.693614e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1.445000e+4 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 1.769194e+2 == Approx( xs[ elas ].value ) );
      REQUIRE( 8.311453e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1.264000e+5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.830077e+1 == Approx( xs[ elas ].value ) );
      REQUIRE( 1.278686e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1.504700e+5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 3.301758e+1 == Approx( xs[ elas ].value ) );
      REQUIRE( 1.061047e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1.779400e+5 * electronVolt, xs );
      REQUIRE( 2 == xs.size() );
      REQUIRE( 2.010266e+1 == Approx( xs[ elas ].value ) );
      REQUIRE( 1.099364e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

    }
  } // GIVEN
} // SCENARIO


