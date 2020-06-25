#include "catch.hpp"
#include "resonanceReconstruction.hpp"

#include <iomanip>

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
    ParticlePair in( neutron, fe54 );
    ParticlePair out( photon, fe55 );

    // channels
    Channel< Photon > capture( in, out, 0. * electronVolt, { 0, 0.0, 0.5, +1 },
                               { 0.0 * rootBarn },
                               0.0 );
    Channel< Neutron > elastic( in, in, 0. * electronVolt, { 0, 0.5, 0.5, +1 },
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
        group1( { elastic }, std::move( single ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group2( { elastic }, std::move( multiple ) );
    SpinGroup< ReichMoore, Constant >
        group3( { elastic }, std::move( single2 ) );
    SpinGroup< ReichMoore, Constant >
        group4( { elastic }, std::move( multiple2 ) );

    CompoundSystem< ReichMoore, ShiftFactor > system1( { group1 } );
    CompoundSystem< ReichMoore, ShiftFactor > system2( { group2 } );
    CompoundSystem< ReichMoore, Constant > system3( { group3 } );
    CompoundSystem< ReichMoore, Constant > system4( { group4 } );

    ReactionID t11 = "n,Fe54_e0{0,1/2,1/2+}->n,Fe54_e0{0,1/2,1/2+}";

    THEN( "T matrix elements can be calculated for a single resonance using "
          "the ShiftFactor boundary condition" ) {

      std::map< ReactionID, std::complex< double > > elements;
      group1.evaluateTMatrix( 1e-5 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 2.7315635336518469E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 2.6262482837542134E-10 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e-4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 8.6379624295084206E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 8.8151189804208580E-10 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e-3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 2.7315638663208645E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 3.2977786512512940E-09 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e-2 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 8.6379732207160647E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.5530453738272414E-08 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e-1 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 2.7315982627877194E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.0013329164342604E-07 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e+0 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 8.6390637708612481E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 8.2704531295578247E-07 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e+1 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 2.7350535821045441E-03 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 7.7363958301328618E-06 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e+2 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 8.7496343285849727E-03 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 7.7390050428265502E-05 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e+3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1308780625816941E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 9.8456461232669743E-04 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e+4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -0.27832633292311088 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  8.4737535716503040E-02 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e+5 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -2.3057837263573785E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  5.3212893363871363E-04 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e+6 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -6.7797364132143496E-03 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  4.5971910221365867E-05 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 7.788000e+3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -1.3081573474631189E-12 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  0.99877608598185263 == Approx( elements[ t11 ].imag() ) );
    } // THEN

    THEN( "T matrix elements can be calculated for multiple resonances using "
          "the ShiftFactor boundary condition" ) {

      std::map< ReactionID, std::complex< double > > elements;
      group2.evaluateTMatrix( 1e-5 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1378337586014968E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 2.7196192971732129E-10 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e-4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 9.9227016858217693E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 9.2734337379377304E-10 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e-3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1378340948540082E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 3.6057602193526968E-09 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e-2 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 9.9227126391879009E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.8134876968625867E-08 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e-1 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1378690529477944E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.2467448407834743E-07 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e+0 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 9.9238214374559137E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.0677325006479699E-06 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e+1 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1413825524341621E-03 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.0131159142876881E-05 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e+2 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 1.0036232268077358E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.0158627872086464E-04 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e+3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.5428380271749681E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.2601846957960908E-03 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e+4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -0.26637088548122501 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  7.6970523479355527E-02 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e+5 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -8.8936876610606230E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  7.9754918283514217E-03 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e+6 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -1.4949017156809678E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  2.2353675680864877E-04 == Approx( elements[ t11 ].imag() ) );

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
      REQUIRE( -1.1241471569546157E-07 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  0.99887879243030220 == Approx( elements[ t11 ].imag() ) );
    } // THEN

    THEN( "T matrix elements can be calculated for a single resonance using "
          "the Constant boundary condition" ) {

      std::map< ReactionID, std::complex< double > > elements;
      group3.evaluateTMatrix( 1e-5 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 2.7315635336518469E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 2.6262482837542134E-10 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e-4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 8.6379624295084206E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 8.8151189804208580E-10 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e-3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 2.7315638663208645E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 3.2977786512512940E-09 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e-2 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 8.6379732207160647E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.5530453738272414E-08 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e-1 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 2.7315982627877194E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.0013329164342604E-07 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e+0 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 8.6390637708612481E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 8.2704531295578247E-07 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e+1 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 2.7350535821045441E-03 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 7.7363958301328618E-06 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e+2 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 8.7496343285849727E-03 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 7.7390050428265502E-05 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e+3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1308780625816941E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 9.8456461232669743E-04 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e+4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -0.27832633292311088 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  8.4737535716503040E-02 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e+5 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -2.3057837263573785E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  5.3212893363871363E-04 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e+6 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -6.7797364132143496E-03 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  4.5971910221365867E-05 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 7.788000e+3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -1.3081573474631189E-12 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  0.99877608598185263 == Approx( elements[ t11 ].imag() ) );
    } // THEN

    THEN( "T matrix elements can be calculated for multiple resonances using "
          "the Constant boundary condition" ) {

      std::map< ReactionID, std::complex< double > > elements;
      group4.evaluateTMatrix( 1e-5 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1378337586014968E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 2.7196192971732129E-10 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e-4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 9.9227016858217693E-06 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 9.2734337379377304E-10 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e-3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1378340948540082E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 3.6057602193526968E-09 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e-2 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 9.9227126391879009E-05 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.8134876968625867E-08 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e-1 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1378690529477944E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.2467448407834743E-07 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e+0 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 9.9238214374559137E-04 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.0677325006479699E-06 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e+1 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.1413825524341621E-03 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.0131159142876881E-05 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e+2 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 1.0036232268077358E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.0158627872086464E-04 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e+3 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( 3.5428380271749681E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE( 1.2601846957960908E-03 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e+4 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -0.26637088548122501 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  7.6970523479355527E-02 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e+5 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -8.8936876610606230E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  7.9754918283514217E-03 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e+6 * electronVolt, elements );
      REQUIRE( 1 == elements.size() );
      REQUIRE( -1.4949017156809678E-02 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  2.2353675680864877E-04 == Approx( elements[ t11 ].imag() ) );

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
      REQUIRE( -1.1241471569546157E-07 == Approx( elements[ t11 ].real() ) );
      REQUIRE(  0.99887879243030220 == Approx( elements[ t11 ].imag() ) );
    } // THEN
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
    ParticlePair out( photon, fe55 );
    ParticlePair in( neutron, fe54 );

    // channels
    Channel< Photon > capture1( in, out, 0. * electronVolt, { 0, 0.0, 0.5, +1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic1( in, in, 0. * electronVolt, { 0, 0.5, 0.5, +1 },
                                 { 5.437300e-1 * rootBarn,
                                   5.437300e-1 * rootBarn } );
    Channel< Photon > capture2( in, out, 0. * electronVolt, { 0, 0.0, 0.5, -1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic2( in, in, 0. * electronVolt, { 1, 0.5, 0.5, -1 },
                                 { 5.437300e-1 * rootBarn,
                                   5.437300e-1 * rootBarn } );
    Channel< Photon > capture3( in, out, 0. * electronVolt, { 0, 0.0, 1.5, -1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic3( in, in, 0. * electronVolt, { 1, 0.5, 1.5, -1 },
                                 { 5.437300e-1 * rootBarn,
                                   5.437300e-1 * rootBarn } );
    Channel< Photon > capture4( in, out, 0. * electronVolt, { 0, 0.0, 1.5, +1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic4( in, in, 0. * electronVolt, { 2, 0.5, 1.5, +1 },
                                 { 5.437300e-1 * rootBarn,
                                   5.437300e-1 * rootBarn } );
    Channel< Photon > capture5( in, out, 0. * electronVolt, { 0, 0.0, 2.5, +1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic5( in, in, 0. * electronVolt, { 2, 0.5, 2.5, +1 },
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
        group1( { elastic1 }, std::move( table1 ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group2( { elastic2 }, std::move( table2 ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group3( { elastic3 }, std::move( table3 ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group4( { elastic4 }, std::move( table4 ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group5( { elastic5 }, std::move( table5 ) );

    CompoundSystem< ReichMoore, ShiftFactor >
        system( { group1, group2, group3, group4, group5 } );

    ReactionID s1t11 = "n,Fe54_e0{0,1/2,1/2+}->n,Fe54_e0{0,1/2,1/2+}";
    ReactionID s2t11 = "n,Fe54_e0{1,1/2,1/2-}->n,Fe54_e0{1,1/2,1/2-}";
    ReactionID s3t11 = "n,Fe54_e0{1,1/2,3/2-}->n,Fe54_e0{1,1/2,3/2-}";
    ReactionID s4t11 = "n,Fe54_e0{2,1/2,3/2+}->n,Fe54_e0{2,1/2,3/2+}";
    ReactionID s5t11 = "n,Fe54_e0{2,1/2,5/2+}->n,Fe54_e0{2,1/2,5/2+}";

    THEN( "T matrix elements can be calculated for multiple resonances using "
          "the ShiftFactor boundary condition" ) {

      std::map< ReactionID, std::complex< double > > elements;
      system.evaluateTMatrix( 1e-5 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( 3.1378337586014968E-06 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE( 2.7196192971732129E-10 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE( 1.6355976862529672E-18 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE( 1.1442049549103583E-23 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( 6.3343763877494508E-20 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE( 4.3753437793456630E-24 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( 6.7456830912506877E-27 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE( 1.0125688683417289E-31 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE( 1.1131371429893516E-30 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE( 4.2880701979043835E-36 == Approx( elements[ s5t11 ].imag() ) );


      system.evaluateTMatrix( 1e-4 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( 9.9227016858217693E-06 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE( 9.2734337379377304E-10 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE( 5.1722140322685651E-17 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE( 3.6182937792296206E-22 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( 2.0031057361716424E-18 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE( 1.3836052630129014E-22 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( 2.1331723121597401E-24 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE( 3.2020239657182831E-29 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE( 3.5200487221536529E-28 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE( 1.3560068609693606E-33 == Approx( elements[ s5t11 ].imag() ) );

      system.evaluateTMatrix( 1e-3 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( 3.1378340948540082E-05 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE( 3.6057602193526968E-09 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE( 1.6355977141078498E-15 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE( 1.1442049956560033E-20 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( 6.3343778474637828E-17 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE( 4.3753463574351308E-21 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( 6.7456837158503738E-22 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE( 1.0125690561088899E-26 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE( 1.1131371505543046E-25 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE( 4.2880702592191839E-31 == Approx( elements[ s5t11 ].imag() ) );

      system.evaluateTMatrix( 1e-2 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( 9.9227126391879009E-05 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE( 1.8134876968625867E-08 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE( 5.1722149131174402E-14 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE( 3.6182950859942320E-19 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( 2.0031103522061930E-15 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE( 1.3836134157102876E-19 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( 2.1331742873190795E-19 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE( 3.2020299034449731E-24 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE( 3.5200489613784904E-23 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE( 1.3560070548638007E-28 == Approx( elements[ s5t11 ].imag() ) );

      system.evaluateTMatrix( 1e-1 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( 3.1378690529477944E-04 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE( 1.2467448407834743E-07 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE( 1.6356004996013098E-12 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE( 1.1442093107580155E-17 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( 6.3345238232562381E-14 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE( 4.3756041822506832E-18 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( 6.7457461764117629E-17 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE( 1.0125878330960735E-21 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE( 1.1131379070501746E-20 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE( 4.2880763907059978E-26 == Approx( elements[ s5t11 ].imag() ) );

      system.evaluateTMatrix( 1e+0 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( 9.9238214374559137E-04 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE( 1.0677325006479699E-06 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE( 5.1723029996526702E-11 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE( 3.6184498194419697E-16 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( 2.0035720935562114E-12 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE( 1.3844291094197795E-16 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( 2.1333718219976592E-14 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE( 3.2026237649595377E-19 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE( 3.5200728840395766E-18 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE( 1.3560264445255817E-23 == Approx( elements[ s5t11 ].imag() ) );

      system.evaluateTMatrix( 1e+1 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( 3.1413825524341621E-03 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE( 1.0131159142876881E-05 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE( 1.6358791010580586E-09 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE( 1.1448815498570552E-14 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( 6.3491651349158321E-11 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE( 4.4015134271892146E-15 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( 6.7519981652668485E-12 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE( 1.0144686561052261E-16 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE( 1.1132135622443639E-15 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE( 4.2886896090906334E-21 == Approx( elements[ s5t11 ].imag() ) );

      system.evaluateTMatrix( 1e+2 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( 1.0036232268077358E-02 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE( 1.0158627872086464E-04 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE( 5.1811281599814464E-08 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE( 3.6581032460483634E-13 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( 2.0511703919988145E-09 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE( 1.4700884584242102E-13 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( 2.1533145481269832E-09 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE( 3.2633269863786763E-14 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE( 3.5224669245007466E-13 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE( 1.3579677016997753E-18 == Approx( elements[ s5t11 ].imag() ) );

      system.evaluateTMatrix( 1e+3 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( 3.5428380271749681E-02 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE( 1.2601846957960908E-03 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE( 1.6642702346261620E-06 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE( 1.4632263688803451E-11 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( 8.4552775690566567E-08 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE( 9.1017936129986887E-12 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( 7.4428141127558590E-07 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE( 1.2884495638645060E-11 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE( 1.1208355764527879E-10 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE( 4.3508285554155527E-16 == Approx( elements[ s5t11 ].imag() ) );

      system.evaluateTMatrix( 1e+4 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( -0.26637088548122501 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE(  7.6970523479355527E-02 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE(  6.2664601273005691E-05 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE(  4.4649716504001658E-09 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE(  1.5295288692520585E-06 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE(  1.7286345586562378E-10 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( -5.1856172798828597E-05 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE(  5.2934868956044800E-07 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE(  3.7810515988498356E-08 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE(  1.5904557618982096E-13 == Approx( elements[ s5t11 ].imag() ) );

      system.evaluateTMatrix( 1e+5 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( -8.8936876610606230E-02 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE(  7.9754918283514217E-03 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE( -1.6837416429181479E-03 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE(  2.8488485136173107E-06 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( -4.7968467173710906E-06 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE(  3.9390097306127509E-11 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( -7.7647765730125733E-03 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE(  6.0309629123224099E-05 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE(  4.2677502728571844E-05 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE(  2.5704845015953520E-09 == Approx( elements[ s5t11 ].imag() ) );

      system.evaluateTMatrix( 1e+6 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( -1.4949017156809678E-02 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE(  2.2353675680864877E-04 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE( -1.2403468646910467E-03 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE(  1.5389580597581292E-06 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( -6.3892797698420987E-06 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE(  4.2732204923486960E-11 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( -0.13633153093569308 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE(  1.8945230557316256E-02 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE( -3.3960659960740171E-04 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE(  1.1553732568883165E-07 == Approx( elements[ s5t11 ].imag() ) );

      system.evaluateTMatrix( 7.788000e+3 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( 1.9602898346729467E-08 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE( 0.99877608598255618 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE( 4.1121898357468075E-05 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE( 2.0272133755200554E-09 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( 3.8221700723879128E-07 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE( 7.9001230125476292E-11 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( 4.5853502872904344E-04 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE( 2.3974737402603093E-07 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE( 1.9910473549888594E-08 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE( 8.1912446437968597E-14 == Approx( elements[ s5t11 ].imag() ) );

      system.evaluateTMatrix( 5.287200e+4 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE(  5.8082713467448460E-09 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE(  0.99900117112954367 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE(  1.1244742492523366E-02 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE(  1.4003412959959896E-04 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( -4.3313296718328622E-06 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE(  5.2078377273691688E-11 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( -3.4271784129272279E-03 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE(  1.1759076070013406E-05 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE(  3.6183715010088403E-06 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE(  3.6500042440603338E-11 == Approx( elements[ s5t11 ].imag() ) );

      system.evaluateTMatrix( 7.190500e+4 * electronVolt, elements );
      REQUIRE( 5 == elements.size() );
      REQUIRE( -1.1241471569546157E-07 == Approx( elements[ s1t11 ].real() ) );
      REQUIRE(  0.99887879243030220 == Approx( elements[ s1t11 ].imag() ) );
      REQUIRE( -2.7421443038996149E-03 == Approx( elements[ s2t11 ].real() ) );
      REQUIRE(  7.5777921859886048E-06 == Approx( elements[ s2t11 ].imag() ) );
      REQUIRE( -4.4916352381910759E-06 == Approx( elements[ s3t11 ].real() ) );
      REQUIRE(  4.3055209147814802E-11 == Approx( elements[ s3t11 ].imag() ) );
      REQUIRE( -5.0393213728232057E-03 == Approx( elements[ s4t11 ].real() ) );
      REQUIRE(  2.5408917902376525E-05 == Approx( elements[ s4t11 ].imag() ) );
      REQUIRE(  1.0099014080820087E-05 == Approx( elements[ s5t11 ].real() ) );
      REQUIRE(  1.8909887190968213E-10 == Approx( elements[ s5t11 ].imag() ) );
    } // THEN
  } // GIVEN*/
} // SCENARIO
