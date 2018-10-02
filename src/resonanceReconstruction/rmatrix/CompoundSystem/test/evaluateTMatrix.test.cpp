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
using ShiftFactor = rmatrix::ShiftFactor;
using Constant = rmatrix::Constant;
using ReichMoore = rmatrix::ReichMoore;
template < typename Formalism, typename Option > using SpinGroup = rmatrix::SpinGroup< Formalism, Option >;
template < typename Formalism, typename Option > using CompoundSystem = rmatrix::CompoundSystem< Formalism, Option >;

constexpr AtomicMass neutronMass = 1.008664 * daltons;

SCENARIO( "evaluateTMatrix" ) {

  GIVEN( "valid data for a CompoundSystem with only one SpinGroup without "
         "missing J values using the Reich Moore formalism" ) {

    // test based on Fe54 ENDF/B-VIII.0 LRF7 resonance evaluation
    // data given in Gamma = 2 gamma^2 P(Er) so conversion is required
    // cross section values extracted from NJOY2016.39

    // because the orbital angular momentum l = 0 for these SpinGroup,
    // SpinGroup< ReichMoore, ShiftFactor > and SpinGroup< ReichMoore, Constant >
    // should give the same results as B = S = 0

    // using SpinGroup< ReichMoore, ShiftFactor > is equivalent to NJOY2016's
    // LRF7 reconstruction

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

    SpinGroup< ReichMoore, ShiftFactor >
        group1( in, { elastic }, std::move( single ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group2( in, { elastic }, std::move( multiple ) );
    SpinGroup< ReichMoore, Constant >
        group3( in, { elastic }, std::move( single2 ) );
    SpinGroup< ReichMoore, Constant >
        group4( in, { elastic }, std::move( multiple2 ) );

    CompoundSystem< ReichMoore, ShiftFactor > system1( { group1 } );
    CompoundSystem< ReichMoore, ShiftFactor > system2( { group2 } );
    CompoundSystem< ReichMoore, Constant > system3( { group3 } );
    CompoundSystem< ReichMoore, Constant > system4( { group4 } );

    ReactionID t11 = "n,Fe54_e0{0,1/2,1/2+}->n,Fe54_e0{0,1/2,1/2+}";

    THEN( "T matrix elements can be calculated for a single resonance using "
          "the ShiftFactor boundary condition" ) {

      tsl::hopscotch_map< ReactionID, std::complex< double > > elements;
      system1.evaluateTMatrix( 1e-5 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 2.7315635336518473E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 2.6262482837542144E-10 == Approx( elements[ t11 ].imag() ) );

      system1.evaluateTMatrix( 1e-4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 8.6379624295084223E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 8.8151189804208590E-10 == Approx( elements[ t11 ].imag() ) );

      system1.evaluateTMatrix( 1e-3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 2.7315638663208651E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 3.2977786512512953E-09 == Approx( elements[ t11 ].imag() ) );

      system1.evaluateTMatrix( 1e-2 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 8.6379732207160647E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.5530453738272408E-08 == Approx( elements[ t11 ].imag() ) );

      system1.evaluateTMatrix( 1e-1 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 2.7315982627877189E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.0013329164342601E-07 == Approx( elements[ t11 ].imag() ) );

      system1.evaluateTMatrix( 1e+0 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 8.6390637708612524E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 8.2704531295578300E-07 == Approx( elements[ t11 ].imag() ) );

      system1.evaluateTMatrix( 1e+1 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 2.7350535821045441E-03 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 7.7363958301328652E-06 == Approx( elements[ t11 ].imag() ) );

      system1.evaluateTMatrix( 1e+2 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 8.7496343285849744E-03 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 7.7390050428265516E-05 == Approx( elements[ t11 ].imag() ) );

      system1.evaluateTMatrix( 1e+3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1308780625816948E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 9.8456461232669787E-04 == Approx( elements[ t11 ].imag() ) );

      system1.evaluateTMatrix( 1e+4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -0.27832633292311082 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 8.4737535716503054E-02 == Approx( elements[ t11 ].imag() ) );

      system1.evaluateTMatrix( 1e+5 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -2.3057837263573792E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 5.3212893363871374E-04 == Approx( elements[ t11 ].imag() ) );

      system1.evaluateTMatrix( 1e+6 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -6.7797364132143505E-03 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 4.5971910221365888E-05 == Approx( elements[ t11 ].imag() ) );

      system1.evaluateTMatrix( 7.788000e+3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -1.3081573474629121E-12 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 0.99877608598185263 == Approx( elements[ t11 ].imag() ) );
    }

    THEN( "T matrix elements can be calculated for multiple resonances using "
          "the ShiftFactor boundary condition" ) {

      tsl::hopscotch_map< ReactionID, std::complex< double > > elements;
      system2.evaluateTMatrix( 1e-5 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1378337586014977E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 2.7196192971732139E-10 == Approx( elements[ t11 ].imag() ) );

      system2.evaluateTMatrix( 1e-4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 9.9227016858217710E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 9.2734337379377335E-010 == Approx( elements[ t11 ].imag() ) );

      system2.evaluateTMatrix( 1e-3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1378340948540075E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 3.6057602193526968E-09 == Approx( elements[ t11 ].imag() ) );

      system2.evaluateTMatrix( 1e-2 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 9.5272244146304985E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.7035286229900731E-08 == Approx( elements[ t11 ].imag() ) );

      system2.evaluateTMatrix( 1e-1 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1378690529477939E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.2467448407834743E-07 == Approx( elements[ t11 ].imag() ) );

      system2.evaluateTMatrix( 1e+0 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 9.9238214374559159E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.0677325006479706E-06 == Approx( elements[ t11 ].imag() ) );

      system2.evaluateTMatrix( 1e+1 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1413825524341629E-03 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.0131159142876882E-05 == Approx( elements[ t11 ].imag() ) );

      system2.evaluateTMatrix( 1e+2 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 1.0036232268077361E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.0158627872086462E-04 == Approx( elements[ t11 ].imag() ) );

      system2.evaluateTMatrix( 1e+3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.5428380271749681E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.2601846957960910E-03 == Approx( elements[ t11 ].imag() ) );

      system2.evaluateTMatrix( 1e+4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -0.26637088548122506 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 7.6970523479355610E-02 == Approx( elements[ t11 ].imag() ) );

      system2.evaluateTMatrix( 1e+5 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -8.8936876610606244E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 7.9754918283514217E-03 == Approx( elements[ t11 ].imag() ) );

      system2.evaluateTMatrix( 1e+6 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -1.4949017156809682E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 2.2353675680864888E-04 == Approx( elements[ t11 ].imag() ) );

      system2.evaluateTMatrix( 7.788000e+3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 1.9602898346729467E-08 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 0.99877608598255618 == Approx( elements[ t11 ].imag() ) );

      system2.evaluateTMatrix( 5.287200e+4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 5.8082713467448460E-09 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 0.99900117112954367 == Approx( elements[ t11 ].imag() ) );

      system2.evaluateTMatrix( 7.190500e+4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -1.1241471569546159E-07 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 0.99887879243030231 == Approx( elements[ t11 ].imag() ) );
    }

    THEN( "T matrix elements can be calculated for a single resonance using "
          "the Constant boundary condition" ) {

      tsl::hopscotch_map< ReactionID, std::complex< double > > elements;
      system3.evaluateTMatrix( 1e-5 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 2.7315635336518473E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 2.6262482837542144E-10 == Approx( elements[ t11 ].imag() ) );

      system3.evaluateTMatrix( 1e-4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 8.6379624295084223E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 8.8151189804208590E-10 == Approx( elements[ t11 ].imag() ) );

      system3.evaluateTMatrix( 1e-3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 2.7315638663208651E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 3.2977786512512953E-09 == Approx( elements[ t11 ].imag() ) );

      system3.evaluateTMatrix( 1e-2 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 8.6379732207160647E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.5530453738272408E-08 == Approx( elements[ t11 ].imag() ) );

      system3.evaluateTMatrix( 1e-1 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 2.7315982627877189E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.0013329164342601E-07 == Approx( elements[ t11 ].imag() ) );

      system3.evaluateTMatrix( 1e+0 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 8.6390637708612524E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 8.2704531295578300E-07 == Approx( elements[ t11 ].imag() ) );

      system3.evaluateTMatrix( 1e+1 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 2.7350535821045441E-03 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 7.7363958301328652E-06 == Approx( elements[ t11 ].imag() ) );

      system3.evaluateTMatrix( 1e+2 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 8.7496343285849744E-03 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 7.7390050428265516E-05 == Approx( elements[ t11 ].imag() ) );

      system3.evaluateTMatrix( 1e+3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1308780625816948E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 9.8456461232669787E-04 == Approx( elements[ t11 ].imag() ) );

      system3.evaluateTMatrix( 1e+4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -0.27832633292311082 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 8.4737535716503054E-02 == Approx( elements[ t11 ].imag() ) );

      system3.evaluateTMatrix( 1e+5 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -2.3057837263573792E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 5.3212893363871374E-04 == Approx( elements[ t11 ].imag() ) );

      system3.evaluateTMatrix( 1e+6 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -6.7797364132143505E-03 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 4.5971910221365888E-05 == Approx( elements[ t11 ].imag() ) );

      system3.evaluateTMatrix( 7.788000e+3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -1.3081573474629121E-12 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 0.99877608598185263 == Approx( elements[ t11 ].imag() ) );
    }

    THEN( "T matrix elements can be calculated for multiple resonances using "
          "the Constant boundary condition" ) {

      tsl::hopscotch_map< ReactionID, std::complex< double > > elements;
      system4.evaluateTMatrix( 1e-5 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1378337586014977E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 2.7196192971732139E-10 == Approx( elements[ t11 ].imag() ) );

      system4.evaluateTMatrix( 1e-4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 9.9227016858217710E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 9.2734337379377335E-010 == Approx( elements[ t11 ].imag() ) );

      system4.evaluateTMatrix( 1e-3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1378340948540075E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 3.6057602193526968E-09 == Approx( elements[ t11 ].imag() ) );

      system4.evaluateTMatrix( 1e-2 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 9.5272244146304985E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.7035286229900731E-08 == Approx( elements[ t11 ].imag() ) );

      system4.evaluateTMatrix( 1e-1 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1378690529477939E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.2467448407834743E-07 == Approx( elements[ t11 ].imag() ) );

      system4.evaluateTMatrix( 1e+0 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 9.9238214374559159E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.0677325006479706E-06 == Approx( elements[ t11 ].imag() ) );

      system4.evaluateTMatrix( 1e+1 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1413825524341629E-03 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.0131159142876882E-05 == Approx( elements[ t11 ].imag() ) );

      system4.evaluateTMatrix( 1e+2 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 1.0036232268077361E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.0158627872086462E-04 == Approx( elements[ t11 ].imag() ) );

      system4.evaluateTMatrix( 1e+3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.5428380271749681E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.2601846957960910E-03 == Approx( elements[ t11 ].imag() ) );

      system4.evaluateTMatrix( 1e+4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -0.26637088548122506 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 7.6970523479355610E-02 == Approx( elements[ t11 ].imag() ) );

      system4.evaluateTMatrix( 1e+5 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -8.8936876610606244E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 7.9754918283514217E-03 == Approx( elements[ t11 ].imag() ) );

      system4.evaluateTMatrix( 1e+6 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -1.4949017156809682E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 2.2353675680864888E-04 == Approx( elements[ t11 ].imag() ) );

      system4.evaluateTMatrix( 7.788000e+3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 1.9602898346729467E-08 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 0.99877608598255618 == Approx( elements[ t11 ].imag() ) );

      system4.evaluateTMatrix( 5.287200e+4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 5.8082713467448460E-09 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 0.99900117112954367 == Approx( elements[ t11 ].imag() ) );

      system4.evaluateTMatrix( 7.190500e+4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -1.1241471569546159E-07 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 0.99887879243030231 == Approx( elements[ t11 ].imag() ) );
    }
  } // GIVEN

  GIVEN( "valid data for a CompoundSystem with five SpinGroup without "
         "missing J values" ) {

    // test based on Fe54 ENDF/B-VIII.0 LRF7 resonance evaluation
    // data given in Gamma = 2 gamma^2 P(Er) so conversion is required
    // cross section values extracted from NJOY2016.39

    // because the oribital angular momentum values for these SpinGroup are
    // 0, 1 and 2, CompoundSystem< ShiftFactor > and CompoundSystem< Constant > will
    // give different results

    // using SpinGroup< ShiftFactor > is equivalent to NJOY2016's LRF7 reconstruction

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

    SpinGroup< ReichMoore, ShiftFactor >
        group1( in, { elastic1 }, std::move( table1 ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group2( in, { elastic2 }, std::move( table2 ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group3( in, { elastic3 }, std::move( table3 ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group4( in, { elastic4 }, std::move( table4 ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group5( in, { elastic5 }, std::move( table5 ) );

    CompoundSystem< ReichMoore, ShiftFactor >
        system( { group1, group2, group3, group4, group5 } );

    ReactionID s1t11 = "n,Fe54_e0{0,1/2,1/2+}->n,Fe54_e0{0,1/2,1/2+}";
    ReactionID s2t11 = "n,Fe54_e0{1,1/2,1/2-}->n,Fe54_e0{1,1/2,1/2-}";
    ReactionID s3t11 = "n,Fe54_e0{1,1/2,3/2-}->n,Fe54_e0{1,1/2,3/2-}";
    ReactionID s4t11 = "n,Fe54_e0{2,1/2,3/2+}->n,Fe54_e0{2,1/2,3/2+}";
    ReactionID s5t11 = "n,Fe54_e0{2,1/2,5/2+}->n,Fe54_e0{2,1/2,5/2+}";

    THEN( "T matrix elements can be calculated for multiple resonances using "
          "the ShiftFactor boundary condition" ) {

      tsl::hopscotch_map< ReactionID, std::complex< double > > elements;
      system.evaluateTMatrix( 1e-5 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( 3.1378337586014977E-06 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE( 2.7196192971732139E-10 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE( 1.6355950942939346E-18 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE( 1.1442031401912170E-23 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( 6.3343749168413473E-20 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE( 4.3753431949582311E-24 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( 6.7456823238281948E-27 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE( 1.0125687531374899E-31 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE( 1.1131354305043496E-30 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE( 4.2880637334944761E-36 == Approx( elements[ s5t11 ].imag() ) );

      system.evaluateTMatrix( 1e-4 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( 9.9227016858217710E-06 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE( 9.2734337379377335E-10 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE( 5.1722058357744195E-17 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE( 3.6182880405838100E-22 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( 2.0031052710296584E-18 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE( 1.3836050782133651E-22 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( 2.1331720694794372E-24 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE( 3.2020236014104841E-29 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE( 3.5200433068005840E-28 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE( 1.3560048167434538E-33 == Approx( elements[ s5t11 ].imag() ) );

      system.evaluateTMatrix( 1e-3 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( 3.1378340948540075E-05 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE( 3.6057602193526968E-09 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE( 1.6355951221488239E-15 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE( 1.1442031809368324E-20 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( 6.3343763765556828E-17 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE( 4.3753457730475648E-21 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( 6.7456829484278797E-22 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE( 1.0125689409046396E-26 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE( 1.1131354380693023E-25 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE( 4.2880637948092275E-31 == Approx( elements[ s5t11 ].imag() ) );

      system.evaluateTMatrix( 1e-2 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( 9.5272244146304985E-05 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE( 1.7035286229900731E-08 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE( 5.1722067166235207E-14 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE( 3.6182893473474415E-19 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( 2.0031098870642221E-15 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE( 1.3836132309103267E-19 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( 2.1331740446387796E-19 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE( 3.2020295391368424E-24 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE( 3.5200435460254291E-23 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE( 1.3560050106377480E-28 == Approx( elements[ s5t11 ].imag() ) );

      system.evaluateTMatrix( 1e-1 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( 3.1378690529477939E-04 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE( 1.2467448407834743E-07 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE( 1.6355979076429895E-12 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE( 1.1442074960351601E-17 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( 6.3345223523485394E-14 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE( 4.3756035978496702E-18 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( 6.7457454089892720E-17 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE( 1.0125877178907628E-21 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE( 1.1131361945651857E-20 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE( 4.2880699262913768E-26 == Approx( elements[ s5t11 ].imag() ) );

      system.evaluateTMatrix( 1e+0 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( 9.9238214374559159E-04 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE( 1.0677325006479706E-06 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE( 5.1722948031810625E-11 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE( 3.6184440806207367E-16 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( 2.0035716284155065E-12 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE( 1.3844289245772673E-16 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( 2.1333715793173588E-14 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE( 3.2026234006178569E-19 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE( 3.5200674686869328E-18 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE( 1.3560244002847740E-23 == Approx( elements[ s5t11 ].imag() ) );

      system.evaluateTMatrix( 1e+1 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( 3.1413825524341629E-03 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE( 1.0131159142876882E-05 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE( 1.6355979076429895E-12 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE( 1.1442074960351601E-17 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( 6.3345223523485394E-14 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE( 4.3756035978496702E-18 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( 6.7457454089892720E-17 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE( 1.0125877178907628E-21 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE( 1.1131361945651857E-20 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE( 4.2880699262913768E-26 == Approx( elements[ s5t11 ].imag() ) );

      system.evaluateTMatrix( 1e+2 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( 1.0036232268077361E-02 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE( 1.0158627872086462E-04 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE( 1.6358765091702951E-09 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE( 1.1448797340032219E-14 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( 6.3491636640481728E-11 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE( 4.4015128414376543E-15 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( 6.7519973978444041E-12 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE( 1.0144685407936371E-16 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE( 1.1132118497607013E-15 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE( 4.2886831442094180E-21 == Approx( elements[ s5t11 ].imag() ) );

      system.evaluateTMatrix( 1e+3 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( 3.5428380271749681E-02 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE( 1.2601846957960910E-03 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE( 5.1811199657406127E-08 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE( 3.6580974134104098E-13 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( 2.0511699269847025E-09 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE( 1.4700882691778281E-13 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( 2.1533143054468470E-09 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE( 3.2633266185460417E-14 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE( 3.5224615091900835E-13 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE( 1.3579656559821107E-18 == Approx( elements[ s5t11 ].imag() ) );

      system.evaluateTMatrix( 1e+4 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( -0.26637088548122506 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE(  7.6970523479355610E-02 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE(  6.2664521516202912E-05 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE(  4.4649609689943662E-09 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE(  1.5295284166410137E-06 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE(  1.7286344079122047E-10 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( -5.1856197063681596E-05 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE(  5.2934869161135448E-07 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE(  3.7810461880346313E-08 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE(  1.5904535158457642E-13 == Approx( elements[ s5t11 ].imag() ) );


      system.evaluateTMatrix( 1e+5 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( -8.8936876610606244E-02 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE(  7.9754918283514217E-03 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE( -1.6837436461608971E-03 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE(  2.8488552760035397E-06 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( -4.7968580856496193E-06 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE(  3.9390244975512434E-11 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( -7.7647841957127059E-03 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE(  6.0309747528649327E-05 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE(  4.2677485827070459E-05 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE(  2.5704827988578055E-09 == Approx( elements[ s5t11 ].imag() ) );

      system.evaluateTMatrix( 1e+6 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( -1.4949017156809682E-02 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE(  2.2353675680864888E-04 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE( -1.2403613977238625E-03 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE(  1.5389941178868449E-06 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( -6.3893622435780201E-06 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE(  4.2733283469076269E-11 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( -0.13633310269494839 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE(  1.8945675998923082E-02 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE( -3.3961017642101870E-04 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE(  1.1553975727848812E-07 == Approx( elements[ s5t11 ].imag() ) );

      system.evaluateTMatrix( 7.788000e+3 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( 1.9602898346729467E-08 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE( 0.99877608598255618 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE( 4.1121843211183448E-05 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE( 2.0272083888020970E-09 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( 3.8221669428971980E-07 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE( 7.9001226524706688E-11 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( 4.5853501574052047E-04 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE( 2.3974736135478750E-07 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE( 1.9910444582306306E-08 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE( 8.1912329273910295E-14 == Approx( elements[ s5t11 ].imag() ) );

      system.evaluateTMatrix( 5.287200e+4 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE(  5.8082713467448460E-09 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE(  0.99900117112954367 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE(  1.1244741626881711E-02 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE(  1.4003410987380283E-04 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( -4.3313345862850243E-06 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE(  5.2078456785758597E-11 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( -3.4271799695886093E-03 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE(  1.1759086746240030E-05 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE(  3.6183680397574522E-06 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE(  3.6499995898614942E-11 == Approx( elements[ s5t11 ].imag() ) );

      system.evaluateTMatrix( 7.190500e+4 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( -1.1241471569546159E-07 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE(  0.99887879243030231 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE( -2.7421456127006932E-03 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE(  7.5777993916659057E-06 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( -4.4916426656618338E-06 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE(  4.3055313280331665E-11 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( -5.0393247251443729E-03 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE(  2.5408951699923551E-05 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE(  1.0099006635710086E-05 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE(  1.8909866100798317E-10 == Approx( elements[ s5t11 ].imag() ) );
    }
  } // GIVEN*/
} // SCENARIO


