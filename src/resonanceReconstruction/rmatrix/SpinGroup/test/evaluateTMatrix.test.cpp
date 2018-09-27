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
template < typename Formalism, typename Option > using SpinGroup = rmatrix::SpinGroup< Formalism, Option >;
using ReactionID = rmatrix::ReactionID;
using ShiftFactor = rmatrix::ShiftFactor;
using Constant = rmatrix::Constant;
using ReichMoore = rmatrix::ReichMoore;

constexpr AtomicMass neutronMass = 1.008664 * daltons;

SCENARIO( "evaluateTMatrix" ) {

  //! @todo add test with more than one entrance channel
  //! @todo add test with a charged particle channel

  GIVEN( "valid data for a SpinGroup with one eliminated capture channel "
         "and one elastic channel using the Reich Moore formalism" ) {

    // test based on Fe54 ENDF/B-VIII.0 LRF7 resonance evaluation
    // data given in Gamma = 2 gamma^2 P(Er) so conversion is required
    // T matrix values extracted from the samm.f90 routines in NJOY2016.43

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

    SpinGroup< ReichMoore, ShiftFactor > group1( in, { elastic },
                                                 std::move( single ) );
    SpinGroup< ReichMoore, ShiftFactor > group2( in, { elastic },
                                                 std::move( multiple ) );
    SpinGroup< ReichMoore, Constant > group3( in, { elastic },
                                              std::move( single2 ) );
    SpinGroup< ReichMoore, Constant > group4( in, { elastic },
                                              std::move( multiple2 ) );

    ReactionID t11 = "n,Fe54_e0{0,1/2,1/2+}->n,Fe54_e0{0,1/2,1/2+}";

    THEN( "T matrix elements can be calculated for a single resonance using "
          "the ShiftFactor boundary condition" ) {

      tsl::hopscotch_map< ReactionID, std::complex< double > > elements;
      group1.evaluateTMatrix( 1e-5 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 2.7315635336518473E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 2.6262482837542144E-10 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e-4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 8.6379624295084223E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 8.8151189804208590E-10 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e-3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 2.7315638663208651E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 3.2977786512512953E-09 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e-2 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 8.6379732207160647E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.5530453738272408E-08 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e-1 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 2.7315982627877189E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.0013329164342601E-07 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e+0 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 8.6390637708612524E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 8.2704531295578300E-07 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e+1 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 2.7350535821045441E-03 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 7.7363958301328652E-06 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e+2 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 8.7496343285849744E-03 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 7.7390050428265516E-05 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e+3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1308780625816948E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 9.8456461232669787E-04 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e+4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -0.27832633292311082 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 8.4737535716503054E-02 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e+5 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -2.3057837263573792E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 5.3212893363871374E-04 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e+6 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -6.7797364132143505E-03 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 4.5971910221365888E-05 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 7.788000e+3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -1.3081573474629121E-12 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 0.99877608598185263 == Approx( elements[ t11 ].imag() ) );
    }

    THEN( "T matrix elements can be calculated for multiple resonances using "
          "the ShiftFactor boundary condition" ) {

      tsl::hopscotch_map< ReactionID, std::complex< double > > elements;
      group2.evaluateTMatrix( 1e-5 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1378337586014977E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 2.7196192971732139E-10 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e-4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 9.9227016858217710E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 9.2734337379377335E-010 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e-3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1378340948540075E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 3.6057602193526968E-09 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e-2 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 9.5272244146304985E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.7035286229900731E-08 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e-1 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1378690529477939E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.2467448407834743E-07 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e+0 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 9.9238214374559159E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.0677325006479706E-06 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e+1 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1413825524341629E-03 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.0131159142876882E-05 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e+2 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 1.0036232268077361E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.0158627872086462E-04 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e+3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.5428380271749681E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.2601846957960910E-03 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e+4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -0.26637088548122506 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 7.6970523479355610E-02 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e+5 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -8.8936876610606244E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 7.9754918283514217E-03 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e+6 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -1.4949017156809682E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 2.2353675680864888E-04 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 7.788000e+3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 1.9602898346729467E-08 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 0.99877608598255618 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 5.287200e+4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 5.8082713467448460E-09 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 0.99900117112954367 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 7.190500e+4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -1.1241471569546159E-07 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 0.99887879243030231 == Approx( elements[ t11 ].imag() ) );
    }

    THEN( "T matrix elements can be calculated for a single resonance using "
          "the Constant boundary condition" ) {

      tsl::hopscotch_map< ReactionID, std::complex< double > > elements;
      group3.evaluateTMatrix( 1e-5 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 2.7315635336518473E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 2.6262482837542144E-10 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e-4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 8.6379624295084223E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 8.8151189804208590E-10 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e-3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 2.7315638663208651E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 3.2977786512512953E-09 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e-2 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 8.6379732207160647E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.5530453738272408E-08 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e-1 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 2.7315982627877189E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.0013329164342601E-07 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e+0 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 8.6390637708612524E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 8.2704531295578300E-07 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e+1 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 2.7350535821045441E-03 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 7.7363958301328652E-06 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e+2 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 8.7496343285849744E-03 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 7.7390050428265516E-05 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e+3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1308780625816948E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 9.8456461232669787E-04 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e+4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -0.27832633292311082 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 8.4737535716503054E-02 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e+5 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -2.3057837263573792E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 5.3212893363871374E-04 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e+6 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -6.7797364132143505E-03 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 4.5971910221365888E-05 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 7.788000e+3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -1.3081573474629121E-12 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 0.99877608598185263 == Approx( elements[ t11 ].imag() ) );
    }

    THEN( "T matrix elements can be calculated for multiple resonances using "
          "the Constant boundary condition" ) {

      tsl::hopscotch_map< ReactionID, std::complex< double > > elements;
      group4.evaluateTMatrix( 1e-5 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1378337586014977E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 2.7196192971732139E-10 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e-4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 9.9227016858217710E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 9.2734337379377335E-010 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e-3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1378340948540075E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 3.6057602193526968E-09 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e-2 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 9.5272244146304985E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.7035286229900731E-08 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e-1 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1378690529477939E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.2467448407834743E-07 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e+0 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 9.9238214374559159E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.0677325006479706E-06 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e+1 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1413825524341629E-03 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.0131159142876882E-05 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e+2 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 1.0036232268077361E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.0158627872086462E-04 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e+3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.5428380271749681E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.2601846957960910E-03 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e+4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -0.26637088548122506 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 7.6970523479355610E-02 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e+5 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -8.8936876610606244E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 7.9754918283514217E-03 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e+6 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -1.4949017156809682E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 2.2353675680864888E-04 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 7.788000e+3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 1.9602898346729467E-08 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 0.99877608598255618 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 5.287200e+4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 5.8082713467448460E-09 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 0.99900117112954367 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 7.190500e+4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -1.1241471569546159E-07 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 0.99887879243030231 == Approx( elements[ t11 ].imag() ) );
    }
  } // GIVEN

  GIVEN( "valid data for a SpinGroup with one eliminated capture channel, "
         "one elastic channel and two fission channels" ) {

    // test based on Pu239 ENDF/B-VIII.0 LRF3 resonance evaluation (some of the
    // widths used in this test were negative in the original evaluation)
    // data given in Gamma = 2 gamma^2 P(Er) so conversion is required
    // T matrix values extracted from the samm.f90 routines in NJOY2016.43
    // (note: LRF7 in NJOY2016 doesn't add potential scattering for missing J
    // values)

    // because the oribital angular momentum l = 0 for these SpinGroup,
    // SpinGroup< ReichMoore, ShiftFactor > and SpinGroup< ReichMoore, Constant >
    // should give the same results

    // using SpinGroup< ReichMoore, ShiftFactor > is equivalent to NJOY2016's
    // LRF7 reconstruction

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
    Channel< Photon > capture( out1, { 0, 0.0, 0.0, +1 },
                               { 0.0 * rootBarn },
                               0.0 );
    Channel< Neutron > elastic( in, { 0, 0.5, 0.0, +1 },
                                { 9.410000e-1 * rootBarn },
                                0.0 );
    Channel< Fission > fission1( out2, "fission1", { 0, 0.0, 0.0, +1 },
                                 { 9.410000e-1 * rootBarn },
                                 0.0 );
    Channel< Fission > fission2( out2, "fission2", { 0, 0.0, 0.0, +1 },
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
      { elastic.channelID(), fission1.channelID(), fission2.channelID() },
      { Resonance( 1.541700e+1 * electronVolt,
                   { eGamma( 2.056203e-3, 1.541700e+1 * electronVolt ),
                     fGamma( 1.093928e-6 ),
                     fGamma( 7.550000e-1 ) },
                   cGamma( 4.054259e-2 ) ) } );
    ResonanceTable single2 = single;

    // multiple resonance table
    ResonanceTable multiple(
      { elastic.channelID(), fission1.channelID(), fission2.channelID() },
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

    SpinGroup< ReichMoore, ShiftFactor >
        group1( in, { elastic, fission1, fission2 }, std::move( single ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group2( in, { elastic, fission1, fission2 }, std::move( multiple ) );
    SpinGroup< ReichMoore, Constant >
        group3( in, { elastic, fission1, fission2 }, std::move( single2 ) );
    SpinGroup< ReichMoore, Constant >
        group4( in, { elastic, fission1, fission2 }, std::move( multiple2 ) );

    ReactionID t11 = "n,Pu239_e0{0,1/2,0+}->n,Pu239_e0{0,1/2,0+}";
    ReactionID t12 = "n,Pu239_e0{0,1/2,0+}->fission1{0,0,0+}";
    ReactionID t13 = "n,Pu239_e0{0,1/2,0+}->fission2{0,0,0+}";
    ReactionID t21 = "fission1{0,0,0+}->n,Pu239_e0{0,1/2,0+}";
    ReactionID t22 = "fission1{0,0,0+}->fission1{0,0,0+}";
    ReactionID t23 = "fission1{0,0,0+}->fission2{0,0,0+}";
    ReactionID t31 = "fission2{0,0,0+}->n,Pu239_e0{0,1/2,0+}";
    ReactionID t32 = "fission2{0,0,0+}->fission1{0,0,0+}";
    ReactionID t33 = "fission2{0,0,0+}->fission2{0,0,0+}";

    THEN( "cross sections can be calculated for a single resonance using the "
          "ShiftFactor boundary condition" ) {

      tsl::hopscotch_map< ReactionID, std::complex< double > > elements;
      group1.evaluateTMatrix( 1e-5 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 5.3671960278418724E-08 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.3847864558992168E-09 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 4.3622324153697221E-08 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 1.1254965041247961E-09 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 3.6239976480108025E-05 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 9.3502507326789021E-07 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 4.3622324153697221E-08 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 1.1254965041247961E-09 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 3.5454400299505883E-08 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 9.1475647772318287E-10 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 2.9454336922612795E-05 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 7.5994926636445372E-07 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 3.6239976480108025E-05 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 9.3502507326789021E-07 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 2.9454336922612795E-05 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 7.5994926636445372E-07 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 2.4469683769066101E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 6.3134058245241256E-04 == Approx( elements[ t33 ].imag() ) );

      group1.evaluateTMatrix( 1e-4 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 1.6972662944694212E-07 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 4.3791500520677800E-09 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 7.7573132643363509E-08 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 2.0014796084808867E-09 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 6.4445179320999697E-05 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 1.6627627102398560E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 7.7573132643363509E-08 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 2.0014796084808867E-09 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 3.5454606785707817E-08 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 9.1477126281006444E-10 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 2.9454508464483777E-05 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 7.5996154932303623E-07 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 6.4445179320999697E-05 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 1.6627627102398560E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 2.9454508464483777E-05 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 7.5996154932303623E-07 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 2.4469826280348803E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 6.3135078672599883E-04 == Approx( elements[ t33 ].imag() ) );

      group1.evaluateTMatrix( 1e-3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 5.3675401113247217E-07 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.3849901146307377E-08 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 1.3795474467599038E-07 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 3.5596558885425543E-09 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 1.1460821493055331E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 2.9572437549072362E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 1.3795474467599038E-07 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 3.5596558885425543E-09 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 3.5456673231866491E-08 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 9.1489100976103598E-10 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 2.9456225199243792E-05 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 7.6006103110835989E-07 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 1.1460821493055331E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 2.9572437549072362E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 2.9456225199243792E-05 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 7.6006103110835989E-07 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 2.4471252486506601E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 6.3143343288551165E-04 == Approx( elements[ t33 ].imag() ) );

      group1.evaluateTMatrix( 1e-2 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 1.6983553119870931E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 4.3850353118746127E-08 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 2.4546518129987473E-07 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 6.3377402846067968E-09 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 2.0392431099384607E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 5.2651839008380049E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 2.4546518129987473E-07 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 6.3377402846067968E-09 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 3.5477355536447543E-08 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 9.1600064898807826E-10 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 2.9473407370210455E-05 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 7.6098288248307639E-07 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 2.0392431099384607E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 5.2651839008380049E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 2.9473407370210455E-05 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 7.6098288248307639E-07 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 2.4485526862844626E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 6.3219927635639551E-04 == Approx( elements[ t33 ].imag() ) );

      group1.evaluateTMatrix( 1e-1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 5.4021847276018600E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.4032018518734498E-07 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 4.3906697669165816E-07 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 1.1404637676351934E-08 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 3.6476224541440958E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 9.4745937813613342E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 4.3906697669165816E-07 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 1.1404637676351934E-08 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 3.5685527197203429E-08 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 9.2692124340637110E-10 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 2.9646349464332679E-05 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 7.7005534921990878E-07 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 3.6476224541440958E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 9.4745937813613342E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 2.9646349464332679E-05 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 7.7005534921990878E-07 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 2.4629201404378166E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 6.3973638006505930E-04 == Approx( elements[ t33 ].imag() ) );

      group1.evaluateTMatrix( 1e+0 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 1.8148061780271335E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 5.0104320258132132E-07 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 8.2945263703656704E-07 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 2.2900054599886333E-08 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 6.8908167184441693E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 1.9024603942904133E-05 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 8.2945263703656704E-07 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 2.2900054599886333E-08 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 3.7909925887227664E-08 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 1.0466412835780584E-09 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 3.1494305935538459E-05 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 8.6951477794573497E-07 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 6.8908167184441693E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 1.9024603942904133E-05 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 3.1494305935538459E-05 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 8.6951477794573497E-07 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 2.6164422196759650E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 7.2236396645979821E-04 == Approx( elements[ t33 ].imag() ) );

      group1.evaluateTMatrix( 1e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 1.5203097691511610E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.1186916195294556E-05 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 3.9074499991002276E-06 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 2.8752242841696715E-07 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 3.2461795379279115E-03 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 2.3886407351019834E-04 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 3.9074499991002276E-06 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 2.8752242841696715E-07 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 1.0042799043508806E-07 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 7.3898065740013609E-09 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 8.3432235258460862E-05 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 6.1392055932365770E-06 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 3.2461795379279115E-03 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 2.3886407351019834E-04 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 8.3432235258460862E-05 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 6.1392055932365770E-06 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 6.9312726960541748E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 5.1002478804555943E-03 == Approx( elements[ t33 ].imag() ) );

      group1.evaluateTMatrix( 1e+2 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -3.0955886124614513E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  1.4653576661553729E-07 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( -4.4740886351904470E-07 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  2.1178977252477513E-09 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -3.7169240762582818E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  1.7594790107887363E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( -4.4740886351904470E-07 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  2.1178977252477513E-09 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( -6.4664500427992813E-09 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  3.0610211269294307E-11 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -5.3721116883905777E-06 == Approx( elements[ t23 ].real() ) );
      REQUIRE(  2.5429945743877634E-08 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -3.7169240762582818E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  1.7594790107887363E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -5.3721116883905777E-06 == Approx( elements[ t32 ].real() ) );
      REQUIRE(  2.5429945743877634E-08 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -4.4629717699094058E-03 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  2.1126353387350176E-05 == Approx( elements[ t33 ].imag() ) );

      group1.evaluateTMatrix( 1e+3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -8.4097609268990106E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  3.4682701528839591E-09 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( -6.8351018912904503E-08 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  2.8188648984846533E-11 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -5.6783753865740234E-05 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  2.3418192314044255E-08 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( -6.8351018912904503E-08 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  2.8188648984846533E-11 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( -5.5552848969690266E-10 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  2.2910554730869460E-13 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -4.6151459811522609E-07 == Approx( elements[ t23 ].real() ) );
      REQUIRE(  1.9033327102599287E-10 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -5.6783753865740234E-05 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  2.3418192314044255E-08 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -4.6151459811522609E-07 == Approx( elements[ t32 ].real() ) );
      REQUIRE(  1.9033327102599287E-10 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -3.8341098291767064E-04 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  1.5812255305473891E-07 == Approx( elements[ t33 ].imag() ) );

      group1.evaluateTMatrix( 2e+3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -5.9004042760189601E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  1.2174383646679523E-09 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( -4.0326043413850601E-08 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  8.3205268742045363E-12 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -3.3501536041607237E-05 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  6.9124170725259833E-09 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( -4.0326043413850601E-08 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  8.3205268742045363E-12 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( -2.7560650107062913E-10 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  5.6866260727081856E-14 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -2.2896471727121253E-07 == Approx( elements[ t23 ].real() ) );
      REQUIRE(  4.7242598629089141E-11 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -3.3501536041607237E-05 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  6.9124170725259833E-09 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -2.2896471727121253E-07 == Approx( elements[ t32 ].real() ) );
      REQUIRE(  4.7242598629089141E-11 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -1.9021627411340149E-04 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  3.9247580141423222E-08 == Approx( elements[ t33 ].imag() ) );

      group1.evaluateTMatrix( 1.541700e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -9.9557762280082612E-15 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 2.5779880786087968E-03 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( -2.2963416433377052E-16 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 5.9462378877927866E-05 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -1.9077242847409717E-13 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 4.9399367268805852E-02 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( -2.2963416433377052E-16 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 5.9462378877927866E-05 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( -5.2966085437830095E-18 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 1.3715247681557698E-06 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -4.4002462678224855E-15 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 1.1394171746724397E-03 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -1.9077242847409717E-13 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 4.9399367268805852E-02 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -4.4002462678224855E-15 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 1.1394171746724397E-03 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -3.6555782926814681E-12 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 0.94658990350152139 == Approx( elements[ t33 ].imag() ) );
    }

    THEN( "cross sections can be calculated for multiple resonances using the "
          "ShiftFactor boundary condition" ) {

      tsl::hopscotch_map< ReactionID, std::complex< double > > elements;
      group2.evaluateTMatrix( 1e-5 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 8.6077679923294326E-08 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.8570343939049700E-09 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 1.2888136111385667E-05 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 9.7894607528028714E-08 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 4.0081962349002548E-05 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 1.1164456176264452E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 1.2888136111385667E-05 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 9.7894607528028714E-08 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 5.9174467940555911E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 3.7061254808378881E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 4.3208176704733495E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 1.4292482833367123E-05 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 4.0081962349002548E-05 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 1.1164456176264452E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 4.3208176704733495E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 1.4292482833367123E-05 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 2.6443852465238640E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 7.3345359723414572E-04 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 1e-4 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 2.7220272617802817E-07 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 5.8725701956782057E-09 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 2.2918751872081955E-05 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 1.7408904822316730E-07 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 7.1277322290083322E-05 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 1.9853875984638998E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 2.2918751872081955E-05 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 1.7408904822316730E-07 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 5.9174580639235266E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 3.7061755644438843E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 4.3208301361143976E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 1.4293712189812666E-05 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 7.1277322290083322E-05 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 1.9853875984638998E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 4.3208301361143976E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 1.4293712189812666E-05 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 2.6444000379628831E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 7.3346529717261484E-04 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 1e-3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 8.6081865742698476E-07 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.8573138778862959E-08 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 4.0756742605136361E-05 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 3.0961750764710787E-07 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 1.2675800823402827E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 3.5310412374359456E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 4.0756742605136361E-05 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 3.0961750764710787E-07 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 5.9175708022794100E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 3.7064308383133275E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 4.3209575043568258E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 1.4298370051136695E-05 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 1.2675800823402827E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 3.5310412374359456E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 4.3209575043568258E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 1.4298370051136695E-05 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 2.6445480939713096E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 7.3355857090884791E-04 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 1e-2 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 2.7233522480381394E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 5.8799731503855405E-08 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 7.2491079343976574E-05 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 5.5098032781778077E-07 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 2.2553600034014534E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 6.2865891242128442E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 7.2491079343976574E-05 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 5.5098032781778077E-07 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 5.9186985408499057E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 3.7082076123237557E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 4.3222401908481596E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 1.4320818280330812E-05 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 2.2553600034014534E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 6.2865891242128442E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 4.3222401908481596E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 1.4320818280330812E-05 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 2.6460300033077304E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 7.3441706653197992E-04 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 1e-1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 8.6503127273711767E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.8794672394205756E-07 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 1.2916257209232690E-04 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 9.8518792030730707E-07 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 4.0330031176166536E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 1.1307201349134586E-05 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 1.2916257209232690E-04 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 9.8518792030730707E-07 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 5.9300000563923538E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 3.7235694751397917E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 4.3351388772216185E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 1.4470029756190969E-05 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 4.0330031176166536E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 1.1307201349134586E-05 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 4.3351388772216185E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 1.4470029756190969E-05 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 2.6609441240777805E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 7.4284391249763656E-04 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 1e+0 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 2.8639910507909064E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 6.6284973829173407E-07 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 2.3429141040536187E-04 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 1.8428769093161262E-06 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 7.5956489245821671E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 2.2581707931633578E-05 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 2.3429141040536187E-04 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 1.8428769093161262E-06 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 6.0454438704396554E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 3.8742163879110527E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 4.4688779003336244E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 1.5816695069991603E-05 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 7.5956489245821671E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 2.2581707931633578E-05 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 4.4688779003336244E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 1.5816695069991603E-05 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 2.8201016434418791E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 8.3491814606862131E-04 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 1e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 1.9431416180186621E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.2666815936613835E-05 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 5.2259368632340708E-04 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 6.5435669882954701E-06 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 3.4199812581358802E-03 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 2.6145017724662367E-04 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 5.2259368632340708E-04 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 6.5435669882954701E-06 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 7.5085348211873047E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 6.0100811547122435E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 6.6260767326451852E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 5.5721401343776834E-05 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 3.4199812581358802E-03 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 2.6145017724662367E-04 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 6.6260767326451852E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 5.5721401343776834E-05 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 7.2131730438568478E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 5.5074649724989677E-03 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 1e+2 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -1.1371918230606258E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  6.7440028937519117E-07 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( -6.3615465189479975E-04 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  3.7343103412946935E-06 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -4.7521348560226043E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  2.8717398712057162E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( -6.3615465189479975E-04 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  3.7343103412946935E-06 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( -5.3265201903090337E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  3.0309992779761852E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -1.9906761702258988E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE(  2.5004165149004856E-06 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -4.7521348560226043E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  2.8717398712057162E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -1.9906761702258988E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE(  2.5004165149004856E-06 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -5.4078818278515389E-03 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  3.0873216863821491E-05 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 1e+3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -2.3356638158054478E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  9.7715058458843184E-09 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( -6.2892983113352560E-05 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  2.1918820914365050E-08 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -6.9651147053680499E-05 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  3.5260535230245288E-08 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( -6.2892983113352560E-05 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  2.1918820914365050E-08 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( -2.9399211525146107E-04 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  9.5134378450556179E-08 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -1.3971263510252607E-05 == Approx( elements[ t23 ].real() ) );
      REQUIRE(  1.5068664218306784E-08 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -6.9651147053680499E-005 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  3.5260535230245288E-08 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -1.3971263510252607E-05 == Approx( elements[ t32 ].real() ) );
      REQUIRE(  1.5068664218306784E-08 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -4.4949738479358695E-04 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  2.1641657866935924E-07 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 2e+3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -1.6226095958907365E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  3.4220039367265144E-09 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( -3.6498856824950382E-05 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  6.3886519761307961E-09 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -4.1026702757550486E-05 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  1.0475988228178090E-08 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( -3.6498856824950382E-05 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  6.3886519761307961E-09 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( -1.4342807016988428E-04 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  2.3034564822960523E-08 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -6.8723140111676468E-06 == Approx( elements[ t23 ].real() ) );
      REQUIRE(  4.0864233816696655E-09 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -4.1026702757550486E-05 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  1.0475988228178090E-08 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -6.8723140111676468E-06 == Approx( elements[ t32 ].real() ) );
      REQUIRE(  4.0864233816696655E-09 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -2.2271662079016282E-04 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  5.3621475639599743E-08 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 1.541700e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE(  4.6855922882881048E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  2.5784243657876950E-03 == Approx( elements[ t11 ].imag() ) );
      REQUIRE(  6.4132524881508788E-04 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  6.5438089890999950E-05 == Approx( elements[ t12 ].imag() ) );
      REQUIRE(  5.9795040686574901E-07 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  4.9399369213564470E-02 == Approx( elements[ t13 ].imag() ) );
      REQUIRE(  6.4132524881508788E-04 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  6.5438089890999950E-05 == Approx( elements[ t21 ].imag() ) );
      REQUIRE(  8.7906104702958187E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  8.3237160796962160E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -2.8267937532975697E-06 == Approx( elements[ t23 ].real() ) );
      REQUIRE(  1.1394282749259659E-03 == Approx( elements[ t23 ].imag() ) );
      REQUIRE(  5.9795040686574901E-07 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  4.9399369213564470E-02 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -2.8267937532975697E-06 == Approx( elements[ t32 ].real() ) );
      REQUIRE(  1.1394282749259659E-03 == Approx( elements[ t32 ].imag() ) );
      REQUIRE(  9.5700608290524391E-06 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  0.94658991700278861 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 3.232700e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE(  1.0268471501034247E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  4.9379169683549627E-03 == Approx( elements[ t11 ].imag() ) );
      REQUIRE(  1.3355065922058540E-03 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  1.2147115890292269E-02 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -2.2388032699726873E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  5.9918063277982860E-02 == Approx( elements[ t13 ].imag() ) );
      REQUIRE(  1.3355065922058540E-03 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  1.2147115890292269E-02 == Approx( elements[ t21 ].imag() ) );
      REQUIRE(  1.6647456754967625E-02 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  3.0081605500153166E-02 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -1.8116766863097097E-03 == Approx( elements[ t23 ].real() ) );
      REQUIRE(  0.14712588608245986 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -2.2388032699726873E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  5.9918063277982860E-02 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -1.8116766863097097E-03 == Approx( elements[ t32 ].real() ) );
      REQUIRE(  0.14712588608245986 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -1.1880436411440707E-03 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  0.72742863159820514 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 4.753400e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -7.6447370161612227E-005 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  8.7754624192493379E-003 == Approx( elements[ t11 ].imag() ) );
      REQUIRE(  4.4048705382630652E-006 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  9.0882947067441849E-002 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -1.1041511969891327E-003 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  6.2176775301615018E-005 == Approx( elements[ t13 ].imag() ) );
      REQUIRE(  4.4048705382630652E-006 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  9.0882947067441849E-002 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( -5.6430729302691567E-007 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  0.94136780546751575 == Approx( elements[ t22 ].imag() ) );
      REQUIRE(  6.5008073784587894E-005 == Approx( elements[ t23 ].real() ) );
      REQUIRE(  4.4999676874436616E-004 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -1.1041511969891327E-003 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  6.2176775301615018E-005 == Approx( elements[ t31 ].imag() ) );
      REQUIRE(  6.5008073784587894E-005 == Approx( elements[ t32 ].real() ) );
      REQUIRE(  4.4999676874436616E-004 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -1.5953769305589553E-002 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  2.6927772660929319E-004 == Approx( elements[ t33 ].imag() ) );
    }

    THEN( "cross sections can be calculated for a single resonance using the "
          "Constant boundary condition" ) {

      tsl::hopscotch_map< ReactionID, std::complex< double > > elements;
      group3.evaluateTMatrix( 1e-5 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 5.3671960278418724E-08 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.3847864558992168E-09 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 4.3622324153697221E-08 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 1.1254965041247961E-09 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 3.6239976480108025E-05 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 9.3502507326789021E-07 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 4.3622324153697221E-08 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 1.1254965041247961E-09 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 3.5454400299505883E-08 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 9.1475647772318287E-10 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 2.9454336922612795E-05 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 7.5994926636445372E-07 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 3.6239976480108025E-05 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 9.3502507326789021E-07 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 2.9454336922612795E-05 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 7.5994926636445372E-07 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 2.4469683769066101E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 6.3134058245241256E-04 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 1e-4 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 1.6972662944694212E-07 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 4.3791500520677800E-09 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 7.7573132643363509E-08 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 2.0014796084808867E-09 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 6.4445179320999697E-05 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 1.6627627102398560E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 7.7573132643363509E-08 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 2.0014796084808867E-09 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 3.5454606785707817E-08 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 9.1477126281006444E-10 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 2.9454508464483777E-05 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 7.5996154932303623E-07 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 6.4445179320999697E-05 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 1.6627627102398560E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 2.9454508464483777E-05 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 7.5996154932303623E-07 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 2.4469826280348803E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 6.3135078672599883E-04 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 1e-3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 5.3675401113247217E-07 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.3849901146307377E-08 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 1.3795474467599038E-07 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 3.5596558885425543E-09 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 1.1460821493055331E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 2.9572437549072362E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 1.3795474467599038E-07 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 3.5596558885425543E-09 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 3.5456673231866491E-08 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 9.1489100976103598E-10 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 2.9456225199243792E-05 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 7.6006103110835989E-07 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 1.1460821493055331E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 2.9572437549072362E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 2.9456225199243792E-05 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 7.6006103110835989E-07 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 2.4471252486506601E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 6.3143343288551165E-04 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 1e-2 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 1.6983553119870931E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 4.3850353118746127E-08 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 2.4546518129987473E-07 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 6.3377402846067968E-09 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 2.0392431099384607E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 5.2651839008380049E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 2.4546518129987473E-07 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 6.3377402846067968E-09 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 3.5477355536447543E-08 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 9.1600064898807826E-10 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 2.9473407370210455E-05 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 7.6098288248307639E-07 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 2.0392431099384607E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 5.2651839008380049E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 2.9473407370210455E-05 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 7.6098288248307639E-07 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 2.4485526862844626E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 6.3219927635639551E-04 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 1e-1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 5.4021847276018600E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.4032018518734498E-07 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 4.3906697669165816E-07 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 1.1404637676351934E-08 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 3.6476224541440958E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 9.4745937813613342E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 4.3906697669165816E-07 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 1.1404637676351934E-08 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 3.5685527197203429E-08 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 9.2692124340637110E-10 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 2.9646349464332679E-05 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 7.7005534921990878E-07 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 3.6476224541440958E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 9.4745937813613342E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 2.9646349464332679E-05 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 7.7005534921990878E-07 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 2.4629201404378166E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 6.3973638006505930E-04 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 1e+0 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 1.8148061780271335E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 5.0104320258132132E-07 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 8.2945263703656704E-07 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 2.2900054599886333E-08 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 6.8908167184441693E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 1.9024603942904133E-05 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 8.2945263703656704E-07 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 2.2900054599886333E-08 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 3.7909925887227664E-08 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 1.0466412835780584E-09 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 3.1494305935538459E-05 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 8.6951477794573497E-07 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 6.8908167184441693E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 1.9024603942904133E-05 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 3.1494305935538459E-05 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 8.6951477794573497E-07 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 2.6164422196759650E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 7.2236396645979821E-04 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 1e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 1.5203097691511610E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.1186916195294556E-05 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 3.9074499991002276E-06 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 2.8752242841696715E-07 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 3.2461795379279115E-03 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 2.3886407351019834E-04 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 3.9074499991002276E-06 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 2.8752242841696715E-07 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 1.0042799043508806E-07 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 7.3898065740013609E-09 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 8.3432235258460862E-05 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 6.1392055932365770E-06 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 3.2461795379279115E-03 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 2.3886407351019834E-04 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 8.3432235258460862E-05 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 6.1392055932365770E-06 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 6.9312726960541748E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 5.1002478804555943E-03 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 1e+2 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -3.0955886124614513E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  1.4653576661553729E-07 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( -4.4740886351904470E-07 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  2.1178977252477513E-09 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -3.7169240762582818E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  1.7594790107887363E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( -4.4740886351904470E-07 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  2.1178977252477513E-09 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( -6.4664500427992813E-09 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  3.0610211269294307E-11 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -5.3721116883905777E-06 == Approx( elements[ t23 ].real() ) );
      REQUIRE(  2.5429945743877634E-08 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -3.7169240762582818E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  1.7594790107887363E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -5.3721116883905777E-06 == Approx( elements[ t32 ].real() ) );
      REQUIRE(  2.5429945743877634E-08 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -4.4629717699094058E-03 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  2.1126353387350176E-05 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 1e+3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -8.4097609268990106E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  3.4682701528839591E-09 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( -6.8351018912904503E-08 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  2.8188648984846533E-11 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -5.6783753865740234E-05 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  2.3418192314044255E-08 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( -6.8351018912904503E-08 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  2.8188648984846533E-11 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( -5.5552848969690266E-10 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  2.2910554730869460E-13 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -4.6151459811522609E-07 == Approx( elements[ t23 ].real() ) );
      REQUIRE(  1.9033327102599287E-10 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -5.6783753865740234E-05 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  2.3418192314044255E-08 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -4.6151459811522609E-07 == Approx( elements[ t32 ].real() ) );
      REQUIRE(  1.9033327102599287E-10 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -3.8341098291767064E-04 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  1.5812255305473891E-07 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 2e+3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -5.9004042760189601E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  1.2174383646679523E-09 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( -4.0326043413850601E-08 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  8.3205268742045363E-12 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -3.3501536041607237E-05 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  6.9124170725259833E-09 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( -4.0326043413850601E-08 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  8.3205268742045363E-12 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( -2.7560650107062913E-10 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  5.6866260727081856E-14 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -2.2896471727121253E-07 == Approx( elements[ t23 ].real() ) );
      REQUIRE(  4.7242598629089141E-11 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -3.3501536041607237E-05 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  6.9124170725259833E-09 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -2.2896471727121253E-07 == Approx( elements[ t32 ].real() ) );
      REQUIRE(  4.7242598629089141E-11 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -1.9021627411340149E-04 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  3.9247580141423222E-08 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 1.541700e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -9.9557762280082612E-15 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 2.5779880786087968E-03 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( -2.2963416433377052E-16 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 5.9462378877927866E-05 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -1.9077242847409717E-13 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 4.9399367268805852E-02 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( -2.2963416433377052E-16 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 5.9462378877927866E-05 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( -5.2966085437830095E-18 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 1.3715247681557698E-06 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -4.4002462678224855E-15 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 1.1394171746724397E-03 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -1.9077242847409717E-13 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 4.9399367268805852E-02 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -4.4002462678224855E-15 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 1.1394171746724397E-03 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -3.6555782926814681E-12 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 0.94658990350152139 == Approx( elements[ t33 ].imag() ) );
    }

    THEN( "cross sections can be calculated for multiple resonances using the "
          "Constant boundary condition" ) {

      tsl::hopscotch_map< ReactionID, std::complex< double > > elements;
      group4.evaluateTMatrix( 1e-5 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 8.6077679923294326E-08 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.8570343939049700E-09 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 1.2888136111385667E-05 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 9.7894607528028714E-08 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 4.0081962349002548E-05 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 1.1164456176264452E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 1.2888136111385667E-05 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 9.7894607528028714E-08 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 5.9174467940555911E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 3.7061254808378881E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 4.3208176704733495E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 1.4292482833367123E-05 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 4.0081962349002548E-05 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 1.1164456176264452E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 4.3208176704733495E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 1.4292482833367123E-05 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 2.6443852465238640E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 7.3345359723414572E-04 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 1e-4 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 2.7220272617802817E-07 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 5.8725701956782057E-09 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 2.2918751872081955E-05 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 1.7408904822316730E-07 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 7.1277322290083322E-05 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 1.9853875984638998E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 2.2918751872081955E-05 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 1.7408904822316730E-07 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 5.9174580639235266E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 3.7061755644438843E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 4.3208301361143976E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 1.4293712189812666E-05 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 7.1277322290083322E-05 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 1.9853875984638998E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 4.3208301361143976E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 1.4293712189812666E-05 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 2.6444000379628831E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 7.3346529717261484E-04 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 1e-3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 8.6081865742698476E-07 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.8573138778862959E-08 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 4.0756742605136361E-05 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 3.0961750764710787E-07 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 1.2675800823402827E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 3.5310412374359456E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 4.0756742605136361E-05 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 3.0961750764710787E-07 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 5.9175708022794100E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 3.7064308383133275E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 4.3209575043568258E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 1.4298370051136695E-05 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 1.2675800823402827E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 3.5310412374359456E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 4.3209575043568258E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 1.4298370051136695E-05 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 2.6445480939713096E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 7.3355857090884791E-04 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 1e-2 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 2.7233522480381394E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 5.8799731503855405E-08 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 7.2491079343976574E-05 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 5.5098032781778077E-07 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 2.2553600034014534E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 6.2865891242128442E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 7.2491079343976574E-05 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 5.5098032781778077E-07 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 5.9186985408499057E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 3.7082076123237557E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 4.3222401908481596E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 1.4320818280330812E-05 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 2.2553600034014534E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 6.2865891242128442E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 4.3222401908481596E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 1.4320818280330812E-05 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 2.6460300033077304E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 7.3441706653197992E-04 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 1e-1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 8.6503127273711767E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.8794672394205756E-07 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 1.2916257209232690E-04 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 9.8518792030730707E-07 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 4.0330031176166536E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 1.1307201349134586E-05 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 1.2916257209232690E-04 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 9.8518792030730707E-07 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 5.9300000563923538E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 3.7235694751397917E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 4.3351388772216185E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 1.4470029756190969E-05 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 4.0330031176166536E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 1.1307201349134586E-05 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 4.3351388772216185E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 1.4470029756190969E-05 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 2.6609441240777805E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 7.4284391249763656E-04 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 1e+0 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 2.8639910507909064E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 6.6284973829173407E-07 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 2.3429141040536187E-04 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 1.8428769093161262E-06 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 7.5956489245821671E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 2.2581707931633578E-05 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 2.3429141040536187E-04 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 1.8428769093161262E-06 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 6.0454438704396554E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 3.8742163879110527E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 4.4688779003336244E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 1.5816695069991603E-05 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 7.5956489245821671E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 2.2581707931633578E-05 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 4.4688779003336244E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 1.5816695069991603E-05 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 2.8201016434418791E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 8.3491814606862131E-04 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 1e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( 1.9431416180186621E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.2666815936613835E-05 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( 5.2259368632340708E-04 == Approx( elements[ t12 ].real() ) );
      REQUIRE( 6.5435669882954701E-06 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( 3.4199812581358802E-03 == Approx( elements[ t13 ].real() ) );
      REQUIRE( 2.6145017724662367E-04 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( 5.2259368632340708E-04 == Approx( elements[ t21 ].real() ) );
      REQUIRE( 6.5435669882954701E-06 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( 7.5085348211873047E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE( 6.0100811547122435E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( 6.6260767326451852E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( 5.5721401343776834E-05 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( 3.4199812581358802E-03 == Approx( elements[ t31 ].real() ) );
      REQUIRE( 2.6145017724662367E-04 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( 6.6260767326451852E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( 5.5721401343776834E-05 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( 7.2131730438568478E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE( 5.5074649724989677E-03 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 1e+2 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -1.1371918230606258E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  6.7440028937519117E-07 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( -6.3615465189479975E-04 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  3.7343103412946935E-06 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -4.7521348560226043E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  2.8717398712057162E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( -6.3615465189479975E-04 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  3.7343103412946935E-06 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( -5.3265201903090337E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  3.0309992779761852E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -1.9906761702258988E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE(  2.5004165149004856E-06 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -4.7521348560226043E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  2.8717398712057162E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -1.9906761702258988E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE(  2.5004165149004856E-06 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -5.4078818278515389E-03 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  3.0873216863821491E-05 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 1e+3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -2.3356638158054478E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  9.7715058458843184E-09 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( -6.2892983113352560E-05 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  2.1918820914365050E-08 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -6.9651147053680499E-05 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  3.5260535230245288E-08 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( -6.2892983113352560E-05 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  2.1918820914365050E-08 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( -2.9399211525146107E-04 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  9.5134378450556179E-08 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -1.3971263510252607E-05 == Approx( elements[ t23 ].real() ) );
      REQUIRE(  1.5068664218306784E-08 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -6.9651147053680499E-005 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  3.5260535230245288E-08 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -1.3971263510252607E-05 == Approx( elements[ t32 ].real() ) );
      REQUIRE(  1.5068664218306784E-08 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -4.4949738479358695E-04 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  2.1641657866935924E-07 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 2e+3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -1.6226095958907365E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  3.4220039367265144E-09 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( -3.6498856824950382E-05 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  6.3886519761307961E-09 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -4.1026702757550486E-05 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  1.0475988228178090E-08 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( -3.6498856824950382E-05 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  6.3886519761307961E-09 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( -1.4342807016988428E-04 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  2.3034564822960523E-08 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -6.8723140111676468E-06 == Approx( elements[ t23 ].real() ) );
      REQUIRE(  4.0864233816696655E-09 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -4.1026702757550486E-05 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  1.0475988228178090E-08 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -6.8723140111676468E-06 == Approx( elements[ t32 ].real() ) );
      REQUIRE(  4.0864233816696655E-09 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -2.2271662079016282E-04 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  5.3621475639599743E-08 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 1.541700e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE(  4.6855922882881048E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  2.5784243657876950E-03 == Approx( elements[ t11 ].imag() ) );
      REQUIRE(  6.4132524881508788E-04 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  6.5438089890999950E-05 == Approx( elements[ t12 ].imag() ) );
      REQUIRE(  5.9795040686574901E-07 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  4.9399369213564470E-02 == Approx( elements[ t13 ].imag() ) );
      REQUIRE(  6.4132524881508788E-04 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  6.5438089890999950E-05 == Approx( elements[ t21 ].imag() ) );
      REQUIRE(  8.7906104702958187E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  8.3237160796962160E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -2.8267937532975697E-06 == Approx( elements[ t23 ].real() ) );
      REQUIRE(  1.1394282749259659E-03 == Approx( elements[ t23 ].imag() ) );
      REQUIRE(  5.9795040686574901E-07 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  4.9399369213564470E-02 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -2.8267937532975697E-06 == Approx( elements[ t32 ].real() ) );
      REQUIRE(  1.1394282749259659E-03 == Approx( elements[ t32 ].imag() ) );
      REQUIRE(  9.5700608290524391E-06 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  0.94658991700278861 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 3.232700e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE(  1.0268471501034247E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  4.9379169683549627E-03 == Approx( elements[ t11 ].imag() ) );
      REQUIRE(  1.3355065922058540E-03 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  1.2147115890292269E-02 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -2.2388032699726873E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  5.9918063277982860E-02 == Approx( elements[ t13 ].imag() ) );
      REQUIRE(  1.3355065922058540E-03 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  1.2147115890292269E-02 == Approx( elements[ t21 ].imag() ) );
      REQUIRE(  1.6647456754967625E-02 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  3.0081605500153166E-02 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -1.8116766863097097E-03 == Approx( elements[ t23 ].real() ) );
      REQUIRE(  0.14712588608245986 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -2.2388032699726873E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  5.9918063277982860E-02 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -1.8116766863097097E-03 == Approx( elements[ t32 ].real() ) );
      REQUIRE(  0.14712588608245986 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -1.1880436411440707E-03 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  0.72742863159820514 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 4.753400e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -7.6447370161612227E-005 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  8.7754624192493379E-003 == Approx( elements[ t11 ].imag() ) );
      REQUIRE(  4.4048705382630652E-006 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  9.0882947067441849E-002 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -1.1041511969891327E-003 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  6.2176775301615018E-005 == Approx( elements[ t13 ].imag() ) );
      REQUIRE(  4.4048705382630652E-006 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  9.0882947067441849E-002 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( -5.6430729302691567E-007 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  0.94136780546751575 == Approx( elements[ t22 ].imag() ) );
      REQUIRE(  6.5008073784587894E-005 == Approx( elements[ t23 ].real() ) );
      REQUIRE(  4.4999676874436616E-004 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -1.1041511969891327E-003 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  6.2176775301615018E-005 == Approx( elements[ t31 ].imag() ) );
      REQUIRE(  6.5008073784587894E-005 == Approx( elements[ t32 ].real() ) );
      REQUIRE(  4.4999676874436616E-004 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -1.5953769305589553E-002 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  2.6927772660929319E-004 == Approx( elements[ t33 ].imag() ) );
    }
  } // GIVEN

  GIVEN( "valid data for a SpinGroup with a resonance at a negative energy" ) {

    // test based on Pu239 ENDF/B-VIII.0 LRF3 resonance evaluation (some of the
    // widths used in this test were negative in the original evaluation)
    // data given in Gamma = 2 gamma^2 P(Er) so conversion is required
    // T matrix values extracted from the samm.f90 routines in NJOY2016.43
    // (note: LRF7 in NJOY2016 doesn't add potential scattering for missing J
    // values)

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
    Channel< Photon > capture( out1, { 0, 0.0, 0.0, +1 },
                               { 0.0 * rootBarn },
                               0.0 );
    Channel< Neutron > elastic( in, { 0, 0.5, 0.0, +1 },
                                { 9.410000e-1 * rootBarn },
                                0.0 );
    Channel< Fission > fission1( out2, "fission1", { 0, 0.0, 0.0, +1 },
                                 { 9.410000e-1 * rootBarn },
                                 0.0 );
    Channel< Fission > fission2( out2, "fission2", { 0, 0.0, 0.0, +1 },
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
      { elastic.channelID(), fission1.channelID(), fission2.channelID() },
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

    SpinGroup< ReichMoore, ShiftFactor >
        group( in, { elastic, fission1, fission2 }, std::move( table ) );

    ReactionID t11 = "n,Pu239_e0{0,1/2,0+}->n,Pu239_e0{0,1/2,0+}";
    ReactionID t12 = "n,Pu239_e0{0,1/2,0+}->fission1{0,0,0+}";
    ReactionID t13 = "n,Pu239_e0{0,1/2,0+}->fission2{0,0,0+}";
    ReactionID t21 = "fission1{0,0,0+}->n,Pu239_e0{0,1/2,0+}";
    ReactionID t22 = "fission1{0,0,0+}->fission1{0,0,0+}";
    ReactionID t23 = "fission1{0,0,0+}->fission2{0,0,0+}";
    ReactionID t31 = "fission2{0,0,0+}->n,Pu239_e0{0,1/2,0+}";
    ReactionID t32 = "fission2{0,0,0+}->fission1{0,0,0+}";
    ReactionID t33 = "fission2{0,0,0+}->fission2{0,0,0+}";

    THEN( "T matrix elements can be calculated" ) {

      tsl::hopscotch_map< ReactionID, std::complex< double > > elements;
      group.evaluateTMatrix( 1e-5 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -2.1263213571644148E-08 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  1.2971372284823746E-09 == Approx( elements[ t11 ].imag() ) );
      REQUIRE(  1.2801222616473468E-05 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  6.7941497903637634E-08 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -3.2396429913820180E-05 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  7.8409475874204870E-07 == Approx( elements[ t13 ].imag() ) );
      REQUIRE(  1.2801222616473468E-05 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  6.7941497903637634E-08 == Approx( elements[ t21 ].imag() ) );
      REQUIRE(  5.9173845275415646E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  3.7013073766943508E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE(  3.7338374143853841E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( -5.8952953266737007E-06 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -3.2396429913820180E-05 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  7.8409475874204870E-07 == Approx( elements[ t31 ].imag() ) );
      REQUIRE(  3.7338374143853841E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( -5.8952953266737007E-06 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -2.2494768086750164E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  5.3990278696026083E-04 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e-4 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -6.7238980019747505E-08 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  4.1018652337798188E-09 == Approx( elements[ t11 ].imag() ) );
      REQUIRE(  2.2764195478870789E-05 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  1.2081854447813629E-07 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -5.7609509202257274E-05 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  1.3943238526046026E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE(  2.2764195478870789E-05 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  1.2081854447813629E-07 == Approx( elements[ t21 ].imag() ) );
      REQUIRE(  5.9173957991666172E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  3.7013569634782119E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE(  3.7338501782878679E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( -5.8961530002853944E-06 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -5.7609509202257274E-05 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  1.3943238526046026E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE(  3.7338501782878679E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( -5.8961530002853944E-06 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -2.2494619781134262E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  5.3989801277359287E-04 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e-3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -2.1259022133504955E-07 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  1.2969815423711923E-08 == Approx( elements[ t11 ].imag() ) );
      REQUIRE(  4.0481898864244864E-05 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  2.1485435148825179E-07 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -1.0243878171068415E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  2.4791858563276496E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE(  4.0481898864244864E-05 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  2.1485435148825179E-07 == Approx( elements[ t21 ].imag() ) );
      REQUIRE(  5.9175085435045054E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  3.7016105628664736E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE(  3.7339789782970444E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( -5.8985972415680078E-06 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -1.0243878171068415E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  2.4791858563276496E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE(  3.7339789782970444E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( -5.8985972415680078E-06 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -2.2493137565948219E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  5.3983476840927208E-04 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e-2 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -6.7106504705850390E-07 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  4.0968590226101194E-08 == Approx( elements[ t11 ].imag() ) );
      REQUIRE(  7.2002338323246500E-05 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  3.8223560913550099E-07 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -1.8203998944443060E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  4.4029696690378751E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE(  7.2002338323246500E-05 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  3.8223560913550099E-07 == Approx( elements[ t21 ].imag() ) );
      REQUIRE(  5.9186363052359274E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  3.7033810048139991E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE(  3.7352708751740388E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( -5.9036440015630698E-06 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -1.8203998944443060E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  4.4029696690378751E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE(  3.7352708751740388E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( -5.9036440015630698E-06 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -2.2478326761901035E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  5.3915393248015267E-04 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e-1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -2.0841956999526716E-006 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  1.2812120650339697E-007 == Approx( elements[ t11 ].imag() ) );
      REQUIRE(  1.2829353810576006E-004 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  6.8306086753656551E-007 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -3.2151211300204829E-004 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  7.7280408405426291E-006 == Approx( elements[ t13 ].imag() ) );
      REQUIRE(  1.2829353810576006E-004 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  6.8306086753656551E-007 == Approx( elements[ t21 ].imag() ) );
      REQUIRE(  5.9299379347602741E-003 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  3.7187122345964958E-005 == Approx( elements[ t22 ].imag() ) );
      REQUIRE(  3.7482247312349001E-004 == Approx( elements[ t23 ].real() ) );
      REQUIRE( -5.8927666460044368E-006 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -3.2151211300204829E-004 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  7.7280408405426291E-006 == Approx( elements[ t31 ].imag() ) );
      REQUIRE(  3.7482247312349001E-004 == Approx( elements[ t32 ].real() ) );
      REQUIRE( -5.8927666460044368E-006 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -2.2331117399884647E-002 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  5.3225711342881816E-004 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e+0 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -5.4462576894535006E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  3.6410337482659467E-07 == Approx( elements[ t11 ].imag() ) );
      REQUIRE(  2.3274167939792254E-04 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  1.2778852474274565E-06 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -5.3469862330085452E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  1.2091349386787629E-05 == Approx( elements[ t13 ].imag() ) );
      REQUIRE(  2.3274167939792254E-04 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  1.2778852474274565E-06 == Approx( elements[ t21 ].imag() ) );
      REQUIRE(  6.0453823311809726E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  3.8691328433044215E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE(  3.8803017114425175E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( -5.5973232137600732E-06 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -5.3469862330085452E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  1.2091349386787629E-05 == Approx( elements[ t31 ].imag() ) );
      REQUIRE(  3.8803017114425175E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( -5.5973232137600732E-06 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -2.0942592513586003E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  4.6898959315351622E-04 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE(  9.8349840297660823E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  5.8666424348490797E-07 == Approx( elements[ t11 ].imag() ) );
      REQUIRE(  5.1802814213470117E-04 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  3.8246542457674050E-06 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -5.1930424677988922E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  7.2326367453940497E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE(  5.1802814213470117E-04 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  3.8246542457674050E-06 == Approx( elements[ t21 ].imag() ) );
      REQUIRE(  7.5084534405060448E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  5.9971563672885378E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE(  5.6506828686182450E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( -2.2413699337732195E-06 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -5.1930424677988922E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  7.2326367453940497E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE(  5.6506828686182450E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( -2.2413699337732195E-06 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -1.1985927576449062E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  1.5880238733464967E-04 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e+2 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -1.0544977684396772E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  5.8448538114817993E-07 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( -6.3603553777656706E-04 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  3.7080441601114395E-06 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -3.7592077890816607E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  1.8299654650210204E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( -6.3603553777656706E-04 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  3.7080441601114395E-06 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( -5.3265185499620166E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  3.0309253929260578E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -1.9763727134302527E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE(  2.1855787825610630E-06 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -3.7592077890816607E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  1.8299654650210204E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -1.9763727134302527E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE(  2.1855787825610630E-06 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -4.2156509939678217E-03 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  1.8818919181566627E-05 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e+3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -2.3101268773318971E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  9.5118011008826419E-09 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( -6.2890907629854781E-05 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  2.1876974201514177E-08 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -6.7926861065473792E-05 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  3.3566173559980889E-08 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( -6.2890907629854781E-05 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  2.1876974201514177E-08 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( -2.9399209838940537E-04 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  9.5133715390776343E-08 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -1.3957249577180668E-05 == Approx( elements[ t23 ].real() ) );
      REQUIRE(  1.4786591112844455E-08 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -6.7926861065473792E-05 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  3.3566173559980889E-08 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -1.3957249577180668E-05 == Approx( elements[ t32 ].real() ) );
      REQUIRE(  1.4786591112844455E-08 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -4.3785479040545712E-04 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  2.0537574395780619E-07 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 2e+3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -1.6135825291599929E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  3.3754144197369966E-09 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( -3.6498239878646618E-05 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  6.3815818260073140E-09 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -4.0514160580584954E-05 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  1.0221438614421103E-08 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( -3.6498239878646618E-05 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  6.3815818260073140E-09 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( -1.4342806595392830E-04 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  2.3034470358314634E-08 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -6.8688110895830866E-06 == Approx( elements[ t23 ].real() ) );
      REQUIRE(  4.0463484131077047E-09 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -4.0514160580584954E-05 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  1.0221438614421103E-08 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -6.8688110895830866E-06 == Approx( elements[ t32 ].real() ) );
      REQUIRE(  4.0463484131077047E-09 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -2.1980648961572106E-04 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  5.2232838285882199E-08 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 3.232700e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE(  1.0277826749300973E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  4.9379012722419128E-03 == Approx( elements[ t11 ].imag() ) );
      REQUIRE(  1.3300936419776699E-03 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  1.2147497328949604E-02 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -2.1398288784814650E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  5.9917205271759158E-02 == Approx( elements[ t13 ].imag() ) );
      REQUIRE(  1.3300936419776699E-03 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  1.2147497328949604E-02 == Approx( elements[ t21 ].imag() ) );
      REQUIRE(  1.6957765703829921E-02 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  3.0089528618649664E-02 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -2.3799477972216744E-03 == Approx( elements[ t23 ].real() ) );
      REQUIRE(  0.14712054109643591 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -2.1398288784814650E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  5.9917205271759158E-02 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -2.3799477972216744E-03 == Approx( elements[ t32 ].real() ) );
      REQUIRE(  0.14712054109643591 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -1.4763497734518113E-04 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  0.72742164346961768 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 4.753400e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -4.9533155268123813E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  8.7747188920056698E-03 == Approx( elements[ t11 ].imag() ) );
      REQUIRE(  1.7814366984730126E-06 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  9.0883002678005143E-02 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -7.1062891323320850E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  5.1380747988416918E-05 == Approx( elements[ t13 ].imag() ) );
      REQUIRE(  1.7814366984730126E-06 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  9.0883002678005143E-02 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( -3.0860148984907700E-07 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  0.94136780169043976 == Approx( elements[ t22 ].imag() ) );
      REQUIRE(  2.6649964127578672E-05 == Approx( elements[ t23 ].real() ) );
      REQUIRE(  4.5080252732947711E-04 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -7.1062891323320850E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  5.1380747988416918E-05 == Approx( elements[ t31 ].imag() ) );
      REQUIRE(  2.6649964127578672E-05 == Approx( elements[ t32 ].real() ) );
      REQUIRE(  4.5080252732947711E-04 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -1.0199940603831531E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  1.1252686951477421E-04 == Approx( elements[ t33 ].imag() ) );
    }
  } // GIVEN

  GIVEN( "valid data for a SpinGroup with resonances using negative widths" ) {

    // test based on Pu239 ENDF/B-VIII.0 LRF3 resonance evaluation (all
    // widths used in this test are from the original evaluation)
    // data given in Gamma = 2 gamma^2 P(Er) so conversion is required
    // T matrix values extracted from the samm.f90 routines in NJOY2016.43
    // (note: LRF7 in NJOY2016 doesn't add potential scattering for missing J
    // values)

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
    Channel< Photon > capture( out1, { 0, 0.0, 0.0, +1 },
                               { 0.0 * rootBarn },
                               0.0 );
    Channel< Neutron > elastic( in, { 0, 0.5, 0.0, +1 },
                                { 9.410000e-1 * rootBarn },
                                0.0 );
    Channel< Fission > fission1( out2, "fission1", { 0, 0.0, 0.0, +1 },
                                 { 9.410000e-1 * rootBarn },
                                 0.0 );
    Channel< Fission > fission2( out2, "fission2", { 0, 0.0, 0.0, +1 },
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
      { elastic.channelID(), fission1.channelID(), fission2.channelID() },
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

    SpinGroup< ReichMoore, ShiftFactor >
        group( in, { elastic, fission1, fission2 }, std::move( table ) );

    ReactionID t11 = "n,Pu239_e0{0,1/2,0+}->n,Pu239_e0{0,1/2,0+}";
    ReactionID t12 = "n,Pu239_e0{0,1/2,0+}->fission1{0,0,0+}";
    ReactionID t13 = "n,Pu239_e0{0,1/2,0+}->fission2{0,0,0+}";
    ReactionID t21 = "fission1{0,0,0+}->n,Pu239_e0{0,1/2,0+}";
    ReactionID t22 = "fission1{0,0,0+}->fission1{0,0,0+}";
    ReactionID t23 = "fission1{0,0,0+}->fission2{0,0,0+}";
    ReactionID t31 = "fission2{0,0,0+}->n,Pu239_e0{0,1/2,0+}";
    ReactionID t32 = "fission2{0,0,0+}->fission1{0,0,0+}";
    ReactionID t33 = "fission2{0,0,0+}->fission2{0,0,0+}";

    THEN( "T matrix elements can be calculated" ) {

      tsl::hopscotch_map< ReactionID, std::complex< double > > elements;
      group.evaluateTMatrix( 1e-5 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE(  8.6094360117330500E-08 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  1.2969579643023378E-09 == Approx( elements[ t11 ].imag() ) );
      REQUIRE(  1.2801915235402997E-05 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  6.5926507156132330E-08 == Approx( elements[ t12 ].imag() ) );
      REQUIRE(  3.2390394989094045E-05 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  8.9681494780835662E-07 == Approx( elements[ t13 ].imag() ) );
      REQUIRE(  1.2801915235402997E-05 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  6.5926507156132330E-08 == Approx( elements[ t21 ].imag() ) );
      REQUIRE(  5.9174467948909766E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  3.7061252598596977E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -4.3208179897630033E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( -1.4291551760903293E-05 == Approx( elements[ t23 ].imag() ) );
      REQUIRE(  3.2390394989094045E-05 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  8.9681494780835662E-07 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -4.3208179897630033E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( -1.4291551760903293E-05 == Approx( elements[ t32 ].imag() ) );
      REQUIRE(  2.6443852496640871E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  7.3345304025029649E-04 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e-4 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE(  2.7225547498043341E-07 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  4.1014378002697183E-09 == Approx( elements[ t11 ].imag() ) );
      REQUIRE(  2.2765426222498603E-05 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  1.1724030555007685E-07 == Approx( elements[ t12 ].imag() ) );
      REQUIRE(  5.7599528678230319E-05 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  1.5948171792593087E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE(  2.2765426222498603E-05 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  1.1724030555007685E-07 == Approx( elements[ t21 ].imag() ) );
      REQUIRE(  5.9174580665652806E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  3.7061748656442602E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -4.3208311458152487E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( -1.4290767857297056E-05 == Approx( elements[ t23 ].imag() ) );
      REQUIRE(  5.7599528678230319E-05 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  1.5948171792593087E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -4.3208311458152487E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( -1.4290767857297056E-05 == Approx( elements[ t32 ].imag() ) );
      REQUIRE(  2.6444000478933142E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  7.3346353581992495E-04 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e-3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE(  8.6098549373132072E-07 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  1.2971844538545411E-08 == Approx( elements[ t11 ].imag() ) );
      REQUIRE(  4.0484071329922845E-05 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  2.0851550150244880E-07 == Approx( elements[ t12 ].imag() ) );
      REQUIRE(  1.0243439585180304E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  2.8364270447639982E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE(  4.0484071329922845E-05 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  2.0851550150244880E-07 == Approx( elements[ t21 ].imag() ) );
      REQUIRE(  5.9175708106342824E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  3.7064286283464467E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -4.3209606977718330E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( -1.4289058530520222E-05 == Approx( elements[ t23 ].imag() ) );
      REQUIRE(  1.0243439585180304E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  2.8364270447639982E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -4.3209606977718330E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( -1.4289058530520222E-05 == Approx( elements[ t32 ].imag() ) );
      REQUIRE(  2.6445481253788445E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  7.3355300054415160E-04 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e-2 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE(  2.7238806360778973E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  4.1071628822804224E-08 == Approx( elements[ t11 ].imag() ) );
      REQUIRE(  7.2005916425102253E-05 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  3.7103339178456461E-07 == Approx( elements[ t12 ].imag() ) );
      REQUIRE(  1.8226982516167097E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  5.0502838646838517E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE(  7.2005916425102253E-05 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  3.7103339178456461E-07 == Approx( elements[ t21 ].imag() ) );
      REQUIRE(  5.9186985672988390E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  3.7082006184398854E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -4.3222503027064743E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( -1.4291349759794567E-05 == Approx( elements[ t23 ].imag() ) );
      REQUIRE(  1.8226982516167097E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  5.0502838646838517E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -4.3222503027064743E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( -1.4291349759794567E-05 == Approx( elements[ t32 ].imag() ) );
      REQUIRE(  2.6460301027705017E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  7.3439943634519788E-04 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e-1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE(  8.6520079376310241E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  1.3139947315384010E-07 == Approx( elements[ t11 ].imag() ) );
      REQUIRE(  1.2829484118662685E-04 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  6.6234773958743636E-07 == Approx( elements[ t12 ].imag() ) );
      REQUIRE(  3.2614691456065846E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  9.0902175262644802E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE(  1.2829484118662685E-04 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  6.6234773958743636E-07 == Approx( elements[ t21 ].imag() ) );
      REQUIRE(  5.9300001409376933E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  3.7235471877830503E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -4.3351712695469921E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( -1.4376112601285967E-05 == Approx( elements[ t23 ].imag() ) );
      REQUIRE(  3.2614691456065846E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  9.0902175262644802E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -4.3351712695469921E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( -1.4376112601285967E-05 == Approx( elements[ t32 ].imag() ) );
      REQUIRE(  2.6609444431195339E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  7.4278767834262715E-04 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e+0 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE(  2.8646113295746248E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  4.6744028655544253E-07 == Approx( elements[ t11 ].imag() ) );
      REQUIRE(  2.3265382264289991E-04 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  1.2141299452405826E-06 == Approx( elements[ t12 ].imag() ) );
      REQUIRE(  6.1843879931493956E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  1.8291330337424252E-05 == Approx( elements[ t13 ].imag() ) );
      REQUIRE(  2.3265382264289991E-04 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  1.2141299452405826E-06 == Approx( elements[ t21 ].imag() ) );
      REQUIRE(  6.0454441690149657E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  3.8741401117764872E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -4.4689947772132762E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( -1.5494916692727437E-05 == Approx( elements[ t23 ].imag() ) );
      REQUIRE(  6.1843879931493956E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  1.8291330337424252E-05 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -4.4689947772132762E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( -1.5494916692727437E-05 == Approx( elements[ t32 ].imag() ) );
      REQUIRE(  2.8201028109277021E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  8.3472384861434264E-04 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE(  1.9449625356215045E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  1.0372444594343622E-05 == Approx( elements[ t11 ].imag() ) );
      REQUIRE(  5.1510765829269119E-04 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  2.1290330431659549E-06 == Approx( elements[ t12 ].imag() ) );
      REQUIRE(  3.0696316664288203E-03 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  2.3494912154491124E-04 == Approx( elements[ t13 ].imag() ) );
      REQUIRE(  5.1510765829269119E-04 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  2.1290330431659549E-06 == Approx( elements[ t21 ].imag() ) );
      REQUIRE(  7.5085394656881415E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  6.0093079962443333E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -6.6289359270729092E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( -5.2355066387279127E-05 == Approx( elements[ t23 ].imag() ) );
      REQUIRE(  3.0696316664288203E-03 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  2.3494912154491124E-04 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -6.6289359270729092E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( -5.2355066387279127E-05 == Approx( elements[ t32 ].imag() ) );
      REQUIRE(  7.2132075903575080E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  5.5052043773717654E-03 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e+2 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -1.1372032929517255E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  5.1934158312836722E-07 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( -6.3526162811888647E-04 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  3.5812399114220902E-06 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -2.6816325859419470E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  1.4113606496224001E-06 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( -6.3526162811888647E-04 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  3.5812399114220902E-06 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( -5.3265203913305274E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  3.0308858408396081E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE(  1.9907302168212498E-04 == Approx( elements[ t23 ].real() ) );
      REQUIRE( -2.0277682885686721E-06 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -2.6816325859419470E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  1.4113606496224001E-06 == Approx( elements[ t31 ].imag() ) );
      REQUIRE(  1.9907302168212498E-04 == Approx( elements[ t32 ].real() ) );
      REQUIRE( -2.0277682885686721E-06 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -5.4078837827793757E-03 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  3.0719306457874704E-05 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e+3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -2.3356639937483160E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  6.8316924521103619E-09 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( -6.2756282536005869E-05 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  2.0285942144168100E-08 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -4.3916352638395594E-05 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  2.0780444778526720E-08 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( -6.2756282536005869E-05 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  2.0285942144168100E-08 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( -2.9399211546199936E-04 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  9.5117202175650674E-08 == Approx( elements[ t22 ].imag() ) );
      REQUIRE(  1.3971269282865116E-05 == Approx( elements[ t23 ].real() ) );
      REQUIRE( -7.9320697666663603E-09 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -4.3916352638395594E-05 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  2.0780444778526720E-08 == Approx( elements[ t31 ].imag() ) );
      REQUIRE(  1.3971269282865116E-05 == Approx( elements[ t32 ].real() ) );
      REQUIRE( -7.9320697666663603E-09 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -4.4949738781199490E-04 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  2.1349394311927587E-07 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 2e+3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -1.6226096271959023E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  2.4077043297474348E-09 == Approx( elements[ t11 ].imag() ) );
      REQUIRE( -3.6418204959627669E-05 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  5.9144855333236295E-09 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -2.5976368082316336E-05 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  6.2188579236122961E-09 == Approx( elements[ t13 ].imag() ) );
      REQUIRE( -3.6418204959627669E-05 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  5.9144855333236295E-09 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( -1.4342807020535596E-04 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  2.3028683929806514E-08 == Approx( elements[ t22 ].imag() ) );
      REQUIRE(  6.8723149961035684E-06 == Approx( elements[ t23 ].real() ) );
      REQUIRE( -1.6429830151187076E-09 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -2.5976368082316336E-05 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  6.2188579236122961E-09 == Approx( elements[ t31 ].imag() ) );
      REQUIRE(  6.8723149961035684E-06 == Approx( elements[ t32 ].real() ) );
      REQUIRE( -1.6429830151187076E-09 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -2.2271662131030396E-04 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  5.2613057057688354E-08 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1.541700e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE(  9.8098208745806375E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  2.5786090151820669E-03 == Approx( elements[ t11 ].imag() ) );
      REQUIRE(  7.1794528627063491E-04 == Approx( elements[ t12 ].real() ) );
      REQUIRE( -5.2664390914965006E-05 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -2.6735454741440945E-05 == Approx( elements[ t13 ].real() ) );
      REQUIRE(  4.9399281630055678E-02 == Approx( elements[ t13 ].imag() ) );
      REQUIRE(  7.1794528627063491E-04 == Approx( elements[ t21 ].real() ) );
      REQUIRE( -5.2664390914965006E-05 == Approx( elements[ t21 ].imag() ) );
      REQUIRE(  8.7907704686553730E-03 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  8.3343078573257701E-05 == Approx( elements[ t22 ].imag() ) );
      REQUIRE( -6.4495186169242103E-05 == Approx( elements[ t23 ].real() ) );
      REQUIRE( -1.1400747932943545E-03 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -2.6735454741440945E-05 == Approx( elements[ t31 ].real() ) );
      REQUIRE(  4.9399281630055678E-02 == Approx( elements[ t31 ].imag() ) );
      REQUIRE( -6.4495186169242103E-05 == Approx( elements[ t32 ].real() ) );
      REQUIRE( -1.1400747932943545E-03 == Approx( elements[ t32 ].imag() ) );
      REQUIRE(  1.2151850651597081E-05 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  0.94658992569665767 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 3.232700e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -2.3172581926204474E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  4.9392463116775323E-03 == Approx( elements[ t11 ].imag() ) );
      REQUIRE(  9.3031529979401332E-04 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  1.2142454746259022E-02 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -5.4616126769468443E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE( -5.9913588141388419E-02 == Approx( elements[ t13 ].imag() ) );
      REQUIRE(  9.3031529979401332E-04 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  1.2142454746259022E-02 == Approx( elements[ t21 ].imag() ) );
      REQUIRE(  1.6657407521629058E-02 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  3.0080787789397490E-02 == Approx( elements[ t22 ].imag() ) );
      REQUIRE(  1.7964258198506292E-03 == Approx( elements[ t23 ].real() ) );
      REQUIRE( -0.14712700890740524 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -5.4616126769468443E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE( -5.9913588141388419E-02 == Approx( elements[ t31 ].imag() ) );
      REQUIRE(  1.7964258198506292E-03 == Approx( elements[ t32 ].real() ) );
      REQUIRE( -0.14712700890740524 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -1.2795632018180745E-03 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  0.72742940272800227 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 4.753400e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
      REQUIRE( -7.6957724678165405E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  8.7745001866921628E-03 == Approx( elements[ t11 ].imag() ) );
      REQUIRE(  3.8148891044013089E-06 == Approx( elements[ t12 ].real() ) );
      REQUIRE(  9.0882977365792025E-02 == Approx( elements[ t12 ].imag() ) );
      REQUIRE( -5.0676531201553218E-04 == Approx( elements[ t13 ].real() ) );
      REQUIRE( -3.5327861901406043E-05 == Approx( elements[ t13 ].imag() ) );
      REQUIRE(  3.8148891044013089E-06 == Approx( elements[ t21 ].real() ) );
      REQUIRE(  9.0882977365792025E-02 == Approx( elements[ t21 ].imag() ) );
      REQUIRE( -4.5182942087408941E-07 == Approx( elements[ t22 ].real() ) );
      REQUIRE(  0.94136780800441544 == Approx( elements[ t22 ].imag() ) );
      REQUIRE(  8.2692928352256167E-05 == Approx( elements[ t23 ].real() ) );
      REQUIRE( -4.5245392366873027E-04 == Approx( elements[ t23 ].imag() ) );
      REQUIRE( -5.0676531201553218E-04 == Approx( elements[ t31 ].real() ) );
      REQUIRE( -3.5327861901406043E-05 == Approx( elements[ t31 ].imag() ) );
      REQUIRE(  8.2692928352256167E-05 == Approx( elements[ t32 ].real() ) );
      REQUIRE( -4.5245392366873027E-04 == Approx( elements[ t32 ].imag() ) );
      REQUIRE( -1.5953943874582939E-02 == Approx( elements[ t33 ].real() ) );
      REQUIRE(  2.6830949543030761E-04 == Approx( elements[ t33 ].imag() ) );
    }
  } // GIVEN
} // SCENARIO


