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
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 1e-4 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 1e-3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 1e-2 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 1e-1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 1e+0 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 1e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 1e+2 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 1e+3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 2e+3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 1.541700e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 3.232700e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 4.753400e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );
    }

    THEN( "cross sections can be calculated for a single resonance using the "
          "ShiftFactor boundary condition" ) {

      tsl::hopscotch_map< ReactionID, std::complex< double > > elements;
      group3.evaluateTMatrix( 1e-5 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 1e-4 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 1e-3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 1e-2 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 1e-1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 1e+0 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 1e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 1e+2 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 1e+3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 2e+3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 1.541700e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );
    }

    THEN( "cross sections can be calculated for multiple resonances using the "
          "ShiftFactor boundary condition" ) {

      tsl::hopscotch_map< ReactionID, std::complex< double > > elements;
      group4.evaluateTMatrix( 1e-5 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 1e-4 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 1e-3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 1e-2 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 1e-1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 1e+0 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 1e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 1e+2 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 1e+3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 2e+3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 1.541700e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 3.232700e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 4.753400e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );
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
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e-4 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e-3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e-2 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e-1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e+0 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e+2 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e+3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 2e+3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1.541700e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 3.232700e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 4.753400e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );
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
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e-4 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e-3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e-2 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e-1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e+0 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e+2 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e+3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 2e+3 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1.541700e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 3.232700e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 4.753400e+1 * electronVolt, elements );
      REQUIRE( 9 == elements.size() );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t11 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t12 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t13 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t21 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t22 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t23 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t31 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t32 ].imag() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].real() ) );
//      REQUIRE( 2.736137e+0 == Approx( elements[ t33 ].imag() ) );
    }
  } // GIVEN
} // SCENARIO


