#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using Particle = rmatrix::Particle;
using ParticlePair = rmatrix::ParticlePair;
using ParticleID = rmatrix::ParticleID;
using ParticlePairID = rmatrix::ParticlePairID;
using Neutron = rmatrix::Neutron;
using Photon = rmatrix::Photon;
using Fission = rmatrix::Fission;
using ChargedParticle = rmatrix::ChargedParticle;
template < typename Type > using Channel = rmatrix::Channel< Type >;
using Resonance = rmatrix::Resonance;
using ResonanceTable = rmatrix::ResonanceTable;
template < typename Formalism, typename Option > using SpinGroup = rmatrix::SpinGroup< Formalism, Option >;
using ReactionChannelID = rmatrix::ReactionChannelID;
using ShiftFactor = rmatrix::ShiftFactor;
using Constant = rmatrix::Constant;
using ReichMoore = rmatrix::ReichMoore;

constexpr AtomicMass neutronMass = 1.008664 * daltons;
constexpr ElectricalCharge elementary = dimwits::constant::elementaryCharge;

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
    Particle photon( ParticleID( "g" ), 0.0 * daltons,
                     0.0 * elementary, 1., +1);
    Particle neutron( ParticleID( "n" ), neutronMass,
                      0.0 * elementary, 0.5, +1);
    Particle fe55( ParticleID( "Fe55_e0" ), 5.446635e+1 * neutronMass,
                   26.0 * elementary, 0.0, +1);
    Particle fe54( ParticleID( "Fe54_e0" ), 5.347624e+1 * neutronMass,
                   26.0 * elementary, 0.0, +1);

    // particle pairs
    ParticlePair in( neutron, fe54 );
    ParticlePair out( photon, fe55 );

    // channels
    Channel< Photon > capture( in, out, 0.0 * electronVolt, { 0, 0.0, 0.5, +1 },
                               { 0.0 * rootBarn },
                               0.0 );
    Channel< Neutron > elastic( in, in, 0.0 * electronVolt, { 0, 0.5, 0.5, +1 },
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

    SpinGroup< ReichMoore, ShiftFactor > group1( { elastic },
                                                 std::move( single ) );
    SpinGroup< ReichMoore, ShiftFactor > group2( { elastic },
                                                 std::move( multiple ) );
    SpinGroup< ReichMoore, Constant > group3( { elastic },
                                              std::move( single2 ) );
    SpinGroup< ReichMoore, Constant > group4( { elastic },
                                              std::move( multiple2 ) );

    ReactionChannelID t11( "n,Fe54{0,1/2,1/2+}->n,Fe54{0,1/2,1/2+}" );

    THEN( "T matrix elements can be calculated for a single resonance using "
          "the ShiftFactor boundary condition" ) {

      std::map< ReactionChannelID, std::complex< double > > elements;
      group1.evaluateTMatrix( 1e-5 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 2.7315635336518469E-06 == Approx( elements[ t11 ].real() ) );
      CHECK( 2.6262482837542134E-10 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e-4 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 8.6379624295084206E-06 == Approx( elements[ t11 ].real() ) );
      CHECK( 8.8151189804208580E-10 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e-3 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 2.7315638663208645E-05 == Approx( elements[ t11 ].real() ) );
      CHECK( 3.2977786512512940E-09 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e-2 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 8.6379732207160647E-05 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.5530453738272414E-08 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e-1 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 2.7315982627877194E-04 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.0013329164342604E-07 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e+0 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 8.6390637708612481E-04 == Approx( elements[ t11 ].real() ) );
      CHECK( 8.2704531295578247E-07 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e+1 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 2.7350535821045441E-03 == Approx( elements[ t11 ].real() ) );
      CHECK( 7.7363958301328618E-06 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e+2 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 8.7496343285849727E-03 == Approx( elements[ t11 ].real() ) );
      CHECK( 7.7390050428265502E-05 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e+3 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 3.1308780625816941E-02 == Approx( elements[ t11 ].real() ) );
      CHECK( 9.8456461232669743E-04 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e+4 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( -0.27832633292311088 == Approx( elements[ t11 ].real() ) );
      CHECK(  8.4737535716503040E-02 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e+5 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( -2.3057837263573785E-02 == Approx( elements[ t11 ].real() ) );
      CHECK(  5.3212893363871363E-04 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 1e+6 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( -6.7797364132143496E-03 == Approx( elements[ t11 ].real() ) );
      CHECK(  4.5971910221365867E-05 == Approx( elements[ t11 ].imag() ) );

      group1.evaluateTMatrix( 7.788000e+3 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( -1.3081573474631189E-12 == Approx( elements[ t11 ].real() ) );
      CHECK(  0.99877608598185263 == Approx( elements[ t11 ].imag() ) );
    } // THEN

    THEN( "T matrix elements can be calculated for multiple resonances using "
          "the ShiftFactor boundary condition" ) {

      std::map< ReactionChannelID, std::complex< double > > elements;
      group2.evaluateTMatrix( 1e-5 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 3.1378337586014968E-06 == Approx( elements[ t11 ].real() ) );
      CHECK( 2.7196192971732129E-10 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e-4 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 9.9227016858217693E-06 == Approx( elements[ t11 ].real() ) );
      CHECK( 9.2734337379377304E-10 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e-3 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 3.1378340948540082E-05 == Approx( elements[ t11 ].real() ) );
      CHECK( 3.6057602193526968E-09 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e-2 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 9.9227126391879009E-05 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.8134876968625867E-08 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e-1 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 3.1378690529477944E-04 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.2467448407834743E-07 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e+0 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 9.9238214374559137E-04 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.0677325006479699E-06 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e+1 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 3.1413825524341621E-03 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.0131159142876881E-05 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e+2 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 1.0036232268077358E-02 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.0158627872086464E-04 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e+3 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 3.5428380271749681E-02 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.2601846957960908E-03 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e+4 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( -0.26637088548122501 == Approx( elements[ t11 ].real() ) );
      CHECK(  7.6970523479355527E-02 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e+5 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( -8.8936876610606230E-02 == Approx( elements[ t11 ].real() ) );
      CHECK(  7.9754918283514217E-03 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 1e+6 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( -1.4949017156809678E-02 == Approx( elements[ t11 ].real() ) );
      CHECK(  2.2353675680864877E-04 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 7.788000e+3 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 1.9602898346729467E-08 == Approx( elements[ t11 ].real() ) );
      CHECK( 0.99877608598255618 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 5.287200e+4 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 5.8082713467448460E-09 == Approx( elements[ t11 ].real() ) );
      CHECK( 0.99900117112954367 == Approx( elements[ t11 ].imag() ) );

      group2.evaluateTMatrix( 7.190500e+4 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( -1.1241471569546157E-07 == Approx( elements[ t11 ].real() ) );
      CHECK(  0.99887879243030220 == Approx( elements[ t11 ].imag() ) );
    } // THEN

    THEN( "T matrix elements can be calculated for a single resonance using "
          "the Constant boundary condition" ) {

      std::map< ReactionChannelID, std::complex< double > > elements;
      group3.evaluateTMatrix( 1e-5 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 2.7315635336518469E-06 == Approx( elements[ t11 ].real() ) );
      CHECK( 2.6262482837542134E-10 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e-4 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 8.6379624295084206E-06 == Approx( elements[ t11 ].real() ) );
      CHECK( 8.8151189804208580E-10 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e-3 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 2.7315638663208645E-05 == Approx( elements[ t11 ].real() ) );
      CHECK( 3.2977786512512940E-09 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e-2 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 8.6379732207160647E-05 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.5530453738272414E-08 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e-1 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 2.7315982627877194E-04 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.0013329164342604E-07 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e+0 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 8.6390637708612481E-04 == Approx( elements[ t11 ].real() ) );
      CHECK( 8.2704531295578247E-07 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e+1 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 2.7350535821045441E-03 == Approx( elements[ t11 ].real() ) );
      CHECK( 7.7363958301328618E-06 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e+2 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 8.7496343285849727E-03 == Approx( elements[ t11 ].real() ) );
      CHECK( 7.7390050428265502E-05 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e+3 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 3.1308780625816941E-02 == Approx( elements[ t11 ].real() ) );
      CHECK( 9.8456461232669743E-04 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e+4 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( -0.27832633292311088 == Approx( elements[ t11 ].real() ) );
      CHECK(  8.4737535716503040E-02 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e+5 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( -2.3057837263573785E-02 == Approx( elements[ t11 ].real() ) );
      CHECK(  5.3212893363871363E-04 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 1e+6 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( -6.7797364132143496E-03 == Approx( elements[ t11 ].real() ) );
      CHECK(  4.5971910221365867E-05 == Approx( elements[ t11 ].imag() ) );

      group3.evaluateTMatrix( 7.788000e+3 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( -1.3081573474631189E-12 == Approx( elements[ t11 ].real() ) );
      CHECK(  0.99877608598185263 == Approx( elements[ t11 ].imag() ) );
    } // THEN

    THEN( "T matrix elements can be calculated for multiple resonances using "
          "the Constant boundary condition" ) {

      std::map< ReactionChannelID, std::complex< double > > elements;
      group4.evaluateTMatrix( 1e-5 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 3.1378337586014968E-06 == Approx( elements[ t11 ].real() ) );
      CHECK( 2.7196192971732129E-10 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e-4 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 9.9227016858217693E-06 == Approx( elements[ t11 ].real() ) );
      CHECK( 9.2734337379377304E-10 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e-3 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 3.1378340948540082E-05 == Approx( elements[ t11 ].real() ) );
      CHECK( 3.6057602193526968E-09 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e-2 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 9.9227126391879009E-05 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.8134876968625867E-08 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e-1 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 3.1378690529477944E-04 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.2467448407834743E-07 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e+0 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 9.9238214374559137E-04 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.0677325006479699E-06 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e+1 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 3.1413825524341621E-03 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.0131159142876881E-05 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e+2 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 1.0036232268077358E-02 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.0158627872086464E-04 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e+3 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 3.5428380271749681E-02 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.2601846957960908E-03 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e+4 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( -0.26637088548122501 == Approx( elements[ t11 ].real() ) );
      CHECK(  7.6970523479355527E-02 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e+5 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( -8.8936876610606230E-02 == Approx( elements[ t11 ].real() ) );
      CHECK(  7.9754918283514217E-03 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 1e+6 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( -1.4949017156809678E-02 == Approx( elements[ t11 ].real() ) );
      CHECK(  2.2353675680864877E-04 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 7.788000e+3 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 1.9602898346729467E-08 == Approx( elements[ t11 ].real() ) );
      CHECK( 0.99877608598255618 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 5.287200e+4 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( 5.8082713467448460E-09 == Approx( elements[ t11 ].real() ) );
      CHECK( 0.99900117112954367 == Approx( elements[ t11 ].imag() ) );

      group4.evaluateTMatrix( 7.190500e+4 * electronVolt, elements );
      CHECK( 1 == elements.size() );
      CHECK( -1.1241471569546157E-07 == Approx( elements[ t11 ].real() ) );
      CHECK(  0.99887879243030220 == Approx( elements[ t11 ].imag() ) );
    } // THEN
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
    Particle photon( ParticleID( "g" ), 0.0 * daltons,
                     0.0 * elementary, 1., +1);
    Particle neutron( ParticleID( "n" ), neutronMass,
                      0.0 * elementary, 0.5, +1);
    Particle pu240( ParticleID( "Pu240_e0" ), 2.379916e+2 * neutronMass,
                    94.0 * elementary, 0.5, +1);
    Particle pu239( ParticleID( "Pu239_e0" ), 2.369986e+2 * neutronMass,
                    94.0 * elementary, 0.5, +1);

    // particle pairs
    ParticlePair in( neutron, pu239 );
    ParticlePair out1( photon, pu240 );
    ParticlePair out2( neutron, pu239, ParticlePairID( "fission" ) );

    // channels
    Channel< Photon > capture( in, out1, 0.0 * electronVolt,
                               { 0, 0.0, 0.0, +1 },
                               { 0.0 * rootBarn },
                               0.0 );
    Channel< Neutron > elastic( in, in, 0.0 * electronVolt,
                                { 0, 0.5, 0.0, +1 },
                                { 9.410000e-1 * rootBarn },
                                0.0 );
    Channel< Fission > fission1( in, out2, "fission1", 0.0 * electronVolt,
                                 { 0, 0.0, 0.0, +1 },
                                 { 9.410000e-1 * rootBarn },
                                 0.0 );
    Channel< Fission > fission2( in, out2, "fission2", 0.0 * electronVolt,
                                 { 0, 0.0, 0.0, +1 },
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
        group1( { elastic, fission1, fission2 }, std::move( single ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group2( { elastic, fission1, fission2 }, std::move( multiple ) );
    SpinGroup< ReichMoore, Constant >
        group3( { elastic, fission1, fission2 }, std::move( single2 ) );
    SpinGroup< ReichMoore, Constant >
        group4( { elastic, fission1, fission2 }, std::move( multiple2 ) );

    ReactionChannelID t11( "n,Pu239{0,1/2,0+}->n,Pu239{0,1/2,0+}" );
    ReactionChannelID t12( "n,Pu239{0,1/2,0+}->fission1{0,0,0+}" );
    ReactionChannelID t13( "n,Pu239{0,1/2,0+}->fission2{0,0,0+}" );
    ReactionChannelID t21( "fission1{0,0,0+}->n,Pu239{0,1/2,0+}" );
    ReactionChannelID t22( "fission1{0,0,0+}->fission1{0,0,0+}" );
    ReactionChannelID t23( "fission1{0,0,0+}->fission2{0,0,0+}" );
    ReactionChannelID t31( "fission2{0,0,0+}->n,Pu239{0,1/2,0+}" );
    ReactionChannelID t32( "fission2{0,0,0+}->fission1{0,0,0+}" );
    ReactionChannelID t33( "fission2{0,0,0+}->fission2{0,0,0+}" );

    THEN( "T matrix elements can be calculated for a single resonance using "
          "the ShiftFactor boundary condition" ) {

      std::map< ReactionChannelID, std::complex< double > > elements;
      group1.evaluateTMatrix( 1e-5 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 5.3671960278418717E-08 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.3847864558992168E-09 == Approx( elements[ t11 ].imag() ) );
      CHECK( 4.3622324153697214E-08 == Approx( elements[ t12 ].real() ) );
      CHECK( 1.1254965041247959E-09 == Approx( elements[ t12 ].imag() ) );
      CHECK( 3.6239976480108025E-05 == Approx( elements[ t13 ].real() ) );
      CHECK( 9.3502507326789042E-07 == Approx( elements[ t13 ].imag() ) );
      CHECK( 4.3622324153697214E-08 == Approx( elements[ t21 ].real() ) );
      CHECK( 1.1254965041247959E-09 == Approx( elements[ t21 ].imag() ) );
      CHECK( 3.5454400299505877E-08 == Approx( elements[ t22 ].real() ) );
      CHECK( 9.1475647772318266E-10 == Approx( elements[ t22 ].imag() ) );
      CHECK( 2.9454336922612795E-05 == Approx( elements[ t23 ].real() ) );
      CHECK( 7.5994926636445372E-07 == Approx( elements[ t23 ].imag() ) );
      CHECK( 3.6239976480108025E-05 == Approx( elements[ t31 ].real() ) );
      CHECK( 9.3502507326789042E-07 == Approx( elements[ t31 ].imag() ) );
      CHECK( 2.9454336922612795E-05 == Approx( elements[ t32 ].real() ) );
      CHECK( 7.5994926636445372E-07 == Approx( elements[ t32 ].imag() ) );
      CHECK( 2.4469683769066101E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 6.3134058245241256E-04 == Approx( elements[ t33 ].imag() ) );

      group1.evaluateTMatrix( 1e-4 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 1.6972662944694210E-07 == Approx( elements[ t11 ].real() ) );
      CHECK( 4.3791500520677809E-09 == Approx( elements[ t11 ].imag() ) );
      CHECK( 7.7573132643363509E-08 == Approx( elements[ t12 ].real() ) );
      CHECK( 2.0014796084808858E-09 == Approx( elements[ t12 ].imag() ) );
      CHECK( 6.4445179320999683E-05 == Approx( elements[ t13 ].real() ) );
      CHECK( 1.6627627102398558E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK( 7.7573132643363509E-08 == Approx( elements[ t21 ].real() ) );
      CHECK( 2.0014796084808858E-09 == Approx( elements[ t21 ].imag() ) );
      CHECK( 3.5454606785707817E-08 == Approx( elements[ t22 ].real() ) );
      CHECK( 9.1477126281006423E-10 == Approx( elements[ t22 ].imag() ) );
      CHECK( 2.9454508464483770E-05 == Approx( elements[ t23 ].real() ) );
      CHECK( 7.5996154932303602E-07 == Approx( elements[ t23 ].imag() ) );
      CHECK( 6.4445179320999683E-05 == Approx( elements[ t31 ].real() ) );
      CHECK( 1.6627627102398558E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK( 2.9454508464483770E-05 == Approx( elements[ t32 ].real() ) );
      CHECK( 7.5996154932303602E-07 == Approx( elements[ t32 ].imag() ) );
      CHECK( 2.4469826280348800E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 6.3135078672599872E-04 == Approx( elements[ t33 ].imag() ) );

      group1.evaluateTMatrix( 1e-3 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 5.3675401113247238E-07 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.3849901146307382E-08 == Approx( elements[ t11 ].imag() ) );
      CHECK( 1.3795474467599038E-07 == Approx( elements[ t12 ].real() ) );
      CHECK( 3.5596558885425563E-09 == Approx( elements[ t12 ].imag() ) );
      CHECK( 1.1460821493055331E-04 == Approx( elements[ t13 ].real() ) );
      CHECK( 2.9572437549072370E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK( 1.3795474467599038E-07 == Approx( elements[ t21 ].real() ) );
      CHECK( 3.5596558885425563E-09 == Approx( elements[ t21 ].imag() ) );
      CHECK( 3.5456673231866484E-08 == Approx( elements[ t22 ].real() ) );
      CHECK( 9.1489100976103619E-10 == Approx( elements[ t22 ].imag() ) );
      CHECK( 2.9456225199243789E-05 == Approx( elements[ t23 ].real() ) );
      CHECK( 7.6006103110835979E-07 == Approx( elements[ t23 ].imag() ) );
      CHECK( 1.1460821493055331E-04 == Approx( elements[ t31 ].real() ) );
      CHECK( 2.9572437549072370E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK( 2.9456225199243789E-05 == Approx( elements[ t32 ].real() ) );
      CHECK( 7.6006103110835979E-07 == Approx( elements[ t32 ].imag() ) );
      CHECK( 2.4471252486506598E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 6.3143343288551165E-04 == Approx( elements[ t33 ].imag() ) );

      group1.evaluateTMatrix( 1e-2 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 1.6983553119870929E-06 == Approx( elements[ t11 ].real() ) );
      CHECK( 4.3850353118746114E-08 == Approx( elements[ t11 ].imag() ) );
      CHECK( 2.4546518129987479E-07 == Approx( elements[ t12 ].real() ) );
      CHECK( 6.3377402846067960E-09 == Approx( elements[ t12 ].imag() ) );
      CHECK( 2.0392431099384604E-04 == Approx( elements[ t13 ].real() ) );
      CHECK( 5.2651839008380032E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK( 2.4546518129987479E-07 == Approx( elements[ t21 ].real() ) );
      CHECK( 6.3377402846067960E-09 == Approx( elements[ t21 ].imag() ) );
      CHECK( 3.5477355536447543E-08 == Approx( elements[ t22 ].real() ) );
      CHECK( 9.1600064898807805E-10 == Approx( elements[ t22 ].imag() ) );
      CHECK( 2.9473407370210455E-05 == Approx( elements[ t23 ].real() ) );
      CHECK( 7.6098288248307618E-07 == Approx( elements[ t23 ].imag() ) );
      CHECK( 2.0392431099384604E-04 == Approx( elements[ t31 ].real() ) );
      CHECK( 5.2651839008380032E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK( 2.9473407370210455E-05 == Approx( elements[ t32 ].real() ) );
      CHECK( 7.6098288248307618E-07 == Approx( elements[ t32 ].imag() ) );
      CHECK( 2.4485526862844626E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 6.3219927635639530E-04 == Approx( elements[ t33 ].imag() ) );

      group1.evaluateTMatrix( 1e-1 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 5.4021847276018575E-06 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.4032018518734493E-07 == Approx( elements[ t11 ].imag() ) );
      CHECK( 4.3906697669165832E-07 == Approx( elements[ t12 ].real() ) );
      CHECK( 1.1404637676351934E-08 == Approx( elements[ t12 ].imag() ) );
      CHECK( 3.6476224541440953E-04 == Approx( elements[ t13 ].real() ) );
      CHECK( 9.4745937813613308E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK( 4.3906697669165832E-07 == Approx( elements[ t21 ].real() ) );
      CHECK( 1.1404637676351934E-08 == Approx( elements[ t21 ].imag() ) );
      CHECK( 3.5685527197203429E-08 == Approx( elements[ t22 ].real() ) );
      CHECK( 9.2692124340637100E-10 == Approx( elements[ t22 ].imag() ) );
      CHECK( 2.9646349464332672E-05 == Approx( elements[ t23 ].real() ) );
      CHECK( 7.7005534921990856E-07 == Approx( elements[ t23 ].imag() ) );
      CHECK( 3.6476224541440953E-04 == Approx( elements[ t31 ].real() ) );
      CHECK( 9.4745937813613308E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK( 2.9646349464332672E-05 == Approx( elements[ t32 ].real() ) );
      CHECK( 7.7005534921990856E-07 == Approx( elements[ t32 ].imag() ) );
      CHECK( 2.4629201404378163E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 6.3973638006505908E-04 == Approx( elements[ t33 ].imag() ) );

      group1.evaluateTMatrix( 1e+0 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 1.8148061780271335E-05 == Approx( elements[ t11 ].real() ) );
      CHECK( 5.0104320258132132E-07 == Approx( elements[ t11 ].imag() ) );
      CHECK( 8.2945263703656714E-07 == Approx( elements[ t12 ].real() ) );
      CHECK( 2.2900054599886333E-08 == Approx( elements[ t12 ].imag() ) );
      CHECK( 6.8908167184441715E-04 == Approx( elements[ t13 ].real() ) );
      CHECK( 1.9024603942904133E-05 == Approx( elements[ t13 ].imag() ) );
      CHECK( 8.2945263703656714E-07 == Approx( elements[ t21 ].real() ) );
      CHECK( 2.2900054599886333E-08 == Approx( elements[ t21 ].imag() ) );
      CHECK( 3.7909925887227664E-08 == Approx( elements[ t22 ].real() ) );
      CHECK( 1.0466412835780584E-09 == Approx( elements[ t22 ].imag() ) );
      CHECK( 3.1494305935538459E-05 == Approx( elements[ t23 ].real() ) );
      CHECK( 8.6951477794573465E-07 == Approx( elements[ t23 ].imag() ) );
      CHECK( 6.8908167184441715E-04 == Approx( elements[ t31 ].real() ) );
      CHECK( 1.9024603942904133E-05 == Approx( elements[ t31 ].imag() ) );
      CHECK( 3.1494305935538459E-05 == Approx( elements[ t32 ].real() ) );
      CHECK( 8.6951477794573465E-07 == Approx( elements[ t32 ].imag() ) );
      CHECK( 2.6164422196759650E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 7.2236396645979800E-04 == Approx( elements[ t33 ].imag() ) );

      group1.evaluateTMatrix( 1e+1 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 1.5203097691511610E-04 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.1186916195294557E-05 == Approx( elements[ t11 ].imag() ) );
      CHECK( 3.9074499991002267E-06 == Approx( elements[ t12 ].real() ) );
      CHECK( 2.8752242841696720E-07 == Approx( elements[ t12 ].imag() ) );
      CHECK( 3.2461795379279106E-03 == Approx( elements[ t13 ].real() ) );
      CHECK( 2.3886407351019826E-04 == Approx( elements[ t13 ].imag() ) );
      CHECK( 3.9074499991002267E-06 == Approx( elements[ t21 ].real() ) );
      CHECK( 2.8752242841696720E-07 == Approx( elements[ t21 ].imag() ) );
      CHECK( 1.0042799043508806E-07 == Approx( elements[ t22 ].real() ) );
      CHECK( 7.3898065740013609E-09 == Approx( elements[ t22 ].imag() ) );
      CHECK( 8.3432235258460835E-05 == Approx( elements[ t23 ].real() ) );
      CHECK( 6.1392055932365745E-06 == Approx( elements[ t23 ].imag() ) );
      CHECK( 3.2461795379279106E-03 == Approx( elements[ t31 ].real() ) );
      CHECK( 2.3886407351019826E-04 == Approx( elements[ t31 ].imag() ) );
      CHECK( 8.3432235258460835E-05 == Approx( elements[ t32 ].real() ) );
      CHECK( 6.1392055932365745E-06 == Approx( elements[ t32 ].imag() ) );
      CHECK( 6.9312726960541721E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 5.1002478804555926E-03 == Approx( elements[ t33 ].imag() ) );

      group1.evaluateTMatrix( 1e+2 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -3.0955886124614519E-05 == Approx( elements[ t11 ].real() ) );
      CHECK(  1.4653576661553732E-07 == Approx( elements[ t11 ].imag() ) );
      CHECK( -4.4740886351904470E-07 == Approx( elements[ t12 ].real() ) );
      CHECK(  2.1178977252477522E-09 == Approx( elements[ t12 ].imag() ) );
      CHECK( -3.7169240762582824E-04 == Approx( elements[ t13 ].real() ) );
      CHECK(  1.7594790107887367E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK( -4.4740886351904470E-07 == Approx( elements[ t21 ].real() ) );
      CHECK(  2.1178977252477522E-09 == Approx( elements[ t21 ].imag() ) );
      CHECK( -6.4664500427992796E-09 == Approx( elements[ t22 ].real() ) );
      CHECK(  3.0610211269294313E-11 == Approx( elements[ t22 ].imag() ) );
      CHECK( -5.3721116883905769E-06 == Approx( elements[ t23 ].real() ) );
      CHECK(  2.5429945743877634E-08 == Approx( elements[ t23 ].imag() ) );
      CHECK( -3.7169240762582824E-04 == Approx( elements[ t31 ].real() ) );
      CHECK(  1.7594790107887367E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK( -5.3721116883905769E-06 == Approx( elements[ t32 ].real() ) );
      CHECK(  2.5429945743877634E-08 == Approx( elements[ t32 ].imag() ) );
      CHECK( -4.4629717699094058E-03 == Approx( elements[ t33 ].real() ) );
      CHECK(  2.1126353387350176E-05 == Approx( elements[ t33 ].imag() ) );

      group1.evaluateTMatrix( 1e+3 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -8.4097609268990072E-06 == Approx( elements[ t11 ].real() ) );
      CHECK(  3.4682701528839583E-09 == Approx( elements[ t11 ].imag() ) );
      CHECK( -6.8351018912904516E-08 == Approx( elements[ t12 ].real() ) );
      CHECK(  2.8188648984846527E-11 == Approx( elements[ t12 ].imag() ) );
      CHECK( -5.6783753865740247E-05 == Approx( elements[ t13 ].real() ) );
      CHECK(  2.3418192314044245E-08 == Approx( elements[ t13 ].imag() ) );
      CHECK( -6.8351018912904516E-08 == Approx( elements[ t21 ].real() ) );
      CHECK(  2.8188648984846527E-11 == Approx( elements[ t21 ].imag() ) );
      CHECK( -5.5552848969690266E-10 == Approx( elements[ t22 ].real() ) );
      CHECK(  2.2910554730869455E-13 == Approx( elements[ t22 ].imag() ) );
      CHECK( -4.6151459811522609E-07 == Approx( elements[ t23 ].real() ) );
      CHECK(  1.9033327102599282E-10 == Approx( elements[ t23 ].imag() ) );
      CHECK( -5.6783753865740247E-05 == Approx( elements[ t31 ].real() ) );
      CHECK(  2.3418192314044245E-08 == Approx( elements[ t31 ].imag() ) );
      CHECK( -4.6151459811522609E-07 == Approx( elements[ t32 ].real() ) );
      CHECK(  1.9033327102599282E-10 == Approx( elements[ t32 ].imag() ) );
      CHECK( -3.8341098291767064E-04 == Approx( elements[ t33 ].real() ) );
      CHECK(  1.5812255305473888E-07 == Approx( elements[ t33 ].imag() ) );

      group1.evaluateTMatrix( 2e+3 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -5.9004042760189593E-06 == Approx( elements[ t11 ].real() ) );
      CHECK(  1.2174383646679529E-09 == Approx( elements[ t11 ].imag() ) );
      CHECK( -4.0326043413850607E-08 == Approx( elements[ t12 ].real() ) );
      CHECK(  8.3205268742045395E-12 == Approx( elements[ t12 ].imag() ) );
      CHECK( -3.3501536041607244E-05 == Approx( elements[ t13 ].real() ) );
      CHECK(  6.9124170725259858E-09 == Approx( elements[ t13 ].imag() ) );
      CHECK( -4.0326043413850607E-08 == Approx( elements[ t21 ].real() ) );
      CHECK(  8.3205268742045395E-12 == Approx( elements[ t21 ].imag() ) );
      CHECK( -2.7560650107062913E-10 == Approx( elements[ t22 ].real() ) );
      CHECK(  5.6866260727081882E-14 == Approx( elements[ t22 ].imag() ) );
      CHECK( -2.2896471727121256E-07 == Approx( elements[ t23 ].real() ) );
      CHECK(  4.7242598629089154E-11 == Approx( elements[ t23 ].imag() ) );
      CHECK( -3.3501536041607244E-05 == Approx( elements[ t31 ].real() ) );
      CHECK(  6.9124170725259858E-09 == Approx( elements[ t31 ].imag() ) );
      CHECK( -2.2896471727121256E-07 == Approx( elements[ t32 ].real() ) );
      CHECK(  4.7242598629089154E-11 == Approx( elements[ t32 ].imag() ) );
      CHECK( -1.9021627411340152E-04 == Approx( elements[ t33 ].real() ) );
      CHECK(  3.9247580141423242E-08 == Approx( elements[ t33 ].imag() ) );

      group1.evaluateTMatrix( 1.541700e+1 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -9.9557762280090690E-15 == Approx( elements[ t11 ].real() ) );
      CHECK(  2.5779880786088033E-03 == Approx( elements[ t11 ].imag() ) );
      CHECK( -2.2963416433378427E-16 == Approx( elements[ t12 ].real() ) );
      CHECK(  5.9462378877928035E-05 == Approx( elements[ t12 ].imag() ) );
      CHECK( -1.9077242847409712E-13 == Approx( elements[ t13 ].real() ) );
      CHECK(  4.9399367268805852E-02 == Approx( elements[ t13 ].imag() ) );
      CHECK( -2.2963416433378427E-16 == Approx( elements[ t21 ].real() ) );
      CHECK(  5.9462378877928035E-05 == Approx( elements[ t21 ].imag() ) );
      CHECK( -5.2966085437830095E-18 == Approx( elements[ t22 ].real() ) );
      CHECK(  1.3715247681557664E-06 == Approx( elements[ t22 ].imag() ) );
      CHECK( -4.4002462678224729E-15 == Approx( elements[ t23 ].real() ) );
      CHECK(  1.1394171746724397E-03 == Approx( elements[ t23 ].imag() ) );
      CHECK( -1.9077242847409712E-13 == Approx( elements[ t31 ].real() ) );
      CHECK(  4.9399367268805852E-02 == Approx( elements[ t31 ].imag() ) );
      CHECK( -4.4002462678224729E-15 == Approx( elements[ t32 ].real() ) );
      CHECK(  1.1394171746724397E-03 == Approx( elements[ t32 ].imag() ) );
      CHECK( -3.6555782926814681E-12 == Approx( elements[ t33 ].real() ) );
      CHECK(  0.94658990350152139 == Approx( elements[ t33 ].imag() ) );
    } // THEN

    THEN( "T matrix elements can be calculated for multiple resonances using "
          "the ShiftFactor boundary condition" ) {

      std::map< ReactionChannelID, std::complex< double > > elements;
      group2.evaluateTMatrix( 1e-5 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 8.6077679923294300E-08 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.8570343939049694E-09 == Approx( elements[ t11 ].imag() ) );
      CHECK( 1.2888136111385671E-05 == Approx( elements[ t12 ].real() ) );
      CHECK( 9.7894607528028767E-08 == Approx( elements[ t12 ].imag() ) );
      CHECK( 4.0081962349002521E-05 == Approx( elements[ t13 ].real() ) );
      CHECK( 1.1164456176264444E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK( 1.2888136111385671E-05 == Approx( elements[ t21 ].real() ) );
      CHECK( 9.7894607528028767E-08 == Approx( elements[ t21 ].imag() ) );
      CHECK( 5.9174467940555911E-03 == Approx( elements[ t22 ].real() ) );
      CHECK( 3.7061254808378901E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK( 4.3208176704733489E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( 1.4292482833367117E-05 == Approx( elements[ t23 ].imag() ) );
      CHECK( 4.0081962349002521E-05 == Approx( elements[ t31 ].real() ) );
      CHECK( 1.1164456176264444E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK( 4.3208176704733489E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( 1.4292482833367117E-05 == Approx( elements[ t32 ].imag() ) );
      CHECK( 2.6443852465238636E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 7.3345359723414551E-04 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 1e-4 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 2.7220272617802817E-07 == Approx( elements[ t11 ].real() ) );
      CHECK( 5.8725701956782081E-09 == Approx( elements[ t11 ].imag() ) );
      CHECK( 2.2918751872081959E-05 == Approx( elements[ t12 ].real() ) );
      CHECK( 1.7408904822316741E-07 == Approx( elements[ t12 ].imag() ) );
      CHECK( 7.1277322290083349E-05 == Approx( elements[ t13 ].real() ) );
      CHECK( 1.9853875984639007E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK( 2.2918751872081959E-05 == Approx( elements[ t21 ].real() ) );
      CHECK( 1.7408904822316741E-07 == Approx( elements[ t21 ].imag() ) );
      CHECK( 5.9174580639235266E-03 == Approx( elements[ t22 ].real() ) );
      CHECK( 3.7061755644438863E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK( 4.3208301361143982E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( 1.4293712189812671E-05 == Approx( elements[ t23 ].imag() ) );
      CHECK( 7.1277322290083349E-05 == Approx( elements[ t31 ].real() ) );
      CHECK( 1.9853875984639007E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK( 4.3208301361143982E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( 1.4293712189812671E-05 == Approx( elements[ t32 ].imag() ) );
      CHECK( 2.6444000379628831E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 7.3346529717261506E-04 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 1e-3 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 8.6081865742698465E-07 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.8573138778862976E-08 == Approx( elements[ t11 ].imag() ) );
      CHECK( 4.0756742605136367E-05 == Approx( elements[ t12 ].real() ) );
      CHECK( 3.0961750764710787E-07 == Approx( elements[ t12 ].imag() ) );
      CHECK( 1.2675800823402835E-04 == Approx( elements[ t13 ].real() ) );
      CHECK( 3.5310412374359472E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK( 4.0756742605136367E-05 == Approx( elements[ t21 ].real() ) );
      CHECK( 3.0961750764710787E-07 == Approx( elements[ t21 ].imag() ) );
      CHECK( 5.9175708022794109E-03 == Approx( elements[ t22 ].real() ) );
      CHECK( 3.7064308383133275E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK( 4.3209575043568269E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( 1.4298370051136698E-05 == Approx( elements[ t23 ].imag() ) );
      CHECK( 1.2675800823402835E-04 == Approx( elements[ t31 ].real() ) );
      CHECK( 3.5310412374359472E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK( 4.3209575043568269E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( 1.4298370051136698E-05 == Approx( elements[ t32 ].imag() ) );
      CHECK( 2.6445480939713099E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 7.3355857090884812E-04 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 1e-2 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 2.7233522480381398E-06 == Approx( elements[ t11 ].real() ) );
      CHECK( 5.8799731503855438E-08 == Approx( elements[ t11 ].imag() ) );
      CHECK( 7.2491079343976574E-05 == Approx( elements[ t12 ].real() ) );
      CHECK( 5.5098032781778119E-07 == Approx( elements[ t12 ].imag() ) );
      CHECK( 2.2553600034014542E-04 == Approx( elements[ t13 ].real() ) );
      CHECK( 6.2865891242128442E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK( 7.2491079343976574E-05 == Approx( elements[ t21 ].real() ) );
      CHECK( 5.5098032781778119E-07 == Approx( elements[ t21 ].imag() ) );
      CHECK( 5.9186985408499057E-03 == Approx( elements[ t22 ].real() ) );
      CHECK( 3.7082076123237578E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK( 4.3222401908481596E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( 1.4320818280330810E-05 == Approx( elements[ t23 ].imag() ) );
      CHECK( 2.2553600034014542E-04 == Approx( elements[ t31 ].real() ) );
      CHECK( 6.2865891242128442E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK( 4.3222401908481596E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( 1.4320818280330810E-05 == Approx( elements[ t32 ].imag() ) );
      CHECK( 2.6460300033077304E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 7.3441706653197981E-04 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 1e-1 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 8.6503127273711784E-06 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.8794672394205753E-07 == Approx( elements[ t11 ].imag() ) );
      CHECK( 1.2916257209232690E-04 == Approx( elements[ t12 ].real() ) );
      CHECK( 9.8518792030730516E-07 == Approx( elements[ t12 ].imag() ) );
      CHECK( 4.0330031176166552E-04 == Approx( elements[ t13 ].real() ) );
      CHECK( 1.1307201349134586E-05 == Approx( elements[ t13 ].imag() ) );
      CHECK( 1.2916257209232690E-04 == Approx( elements[ t21 ].real() ) );
      CHECK( 9.8518792030730516E-07 == Approx( elements[ t21 ].imag() ) );
      CHECK( 5.9300000563923538E-03 == Approx( elements[ t22 ].real() ) );
      CHECK( 3.7235694751397835E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK( 4.3351388772216185E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( 1.4470029756190969E-05 == Approx( elements[ t23 ].imag() ) );
      CHECK( 4.0330031176166552E-04 == Approx( elements[ t31 ].real() ) );
      CHECK( 1.1307201349134586E-05 == Approx( elements[ t31 ].imag() ) );
      CHECK( 4.3351388772216185E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( 1.4470029756190969E-05 == Approx( elements[ t32 ].imag() ) );
      CHECK( 2.6609441240777808E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 7.4284391249763656E-04 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 1e+0 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 2.8639910507909060E-05 == Approx( elements[ t11 ].real() ) );
      CHECK( 6.6284973829173364E-07 == Approx( elements[ t11 ].imag() ) );
      CHECK( 2.3429141040536179E-04 == Approx( elements[ t12 ].real() ) );
      CHECK( 1.8428769093161273E-06 == Approx( elements[ t12 ].imag() ) );
      CHECK( 7.5956489245821639E-04 == Approx( elements[ t13 ].real() ) );
      CHECK( 2.2581707931633565E-05 == Approx( elements[ t13 ].imag() ) );
      CHECK( 2.3429141040536179E-04 == Approx( elements[ t21 ].real() ) );
      CHECK( 1.8428769093161273E-06 == Approx( elements[ t21 ].imag() ) );
      CHECK( 6.0454438704396537E-03 == Approx( elements[ t22 ].real() ) );
      CHECK( 3.8742163879110548E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK( 4.4688779003336233E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( 1.5816695069991596E-05 == Approx( elements[ t23 ].imag() ) );
      CHECK( 7.5956489245821639E-04 == Approx( elements[ t31 ].real() ) );
      CHECK( 2.2581707931633565E-05 == Approx( elements[ t31 ].imag() ) );
      CHECK( 4.4688779003336233E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( 1.5816695069991596E-05 == Approx( elements[ t32 ].imag() ) );
      CHECK( 2.8201016434418785E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 8.3491814606862077E-04 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 1e+1 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 1.9431416180186621E-04 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.2666815936613832E-05 == Approx( elements[ t11 ].imag() ) );
      CHECK( 5.2259368632340708E-04 == Approx( elements[ t12 ].real() ) );
      CHECK( 6.5435669882954769E-06 == Approx( elements[ t12 ].imag() ) );
      CHECK( 3.4199812581358797E-03 == Approx( elements[ t13 ].real() ) );
      CHECK( 2.6145017724662367E-04 == Approx( elements[ t13 ].imag() ) );
      CHECK( 5.2259368632340708E-04 == Approx( elements[ t21 ].real() ) );
      CHECK( 6.5435669882954769E-06 == Approx( elements[ t21 ].imag() ) );
      CHECK( 7.5085348211873038E-03 == Approx( elements[ t22 ].real() ) );
      CHECK( 6.0100811547122537E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK( 6.6260767326451852E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( 5.5721401343776834E-05 == Approx( elements[ t23 ].imag() ) );
      CHECK( 3.4199812581358797E-03 == Approx( elements[ t31 ].real() ) );
      CHECK( 2.6145017724662367E-04 == Approx( elements[ t31 ].imag() ) );
      CHECK( 6.6260767326451852E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( 5.5721401343776834E-05 == Approx( elements[ t32 ].imag() ) );
      CHECK( 7.2131730438568478E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 5.5074649724989677E-03 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 1e+2 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -1.1371918230606257E-04 == Approx( elements[ t11 ].real() ) );
      CHECK(  6.7440028937519096E-07 == Approx( elements[ t11 ].imag() ) );
      CHECK( -6.3615465189479985E-04 == Approx( elements[ t12 ].real() ) );
      CHECK(  3.7343103412946918E-06 == Approx( elements[ t12 ].imag() ) );
      CHECK( -4.7521348560226043E-04 == Approx( elements[ t13 ].real() ) );
      CHECK(  2.8717398712057162E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK( -6.3615465189479985E-04 == Approx( elements[ t21 ].real() ) );
      CHECK(  3.7343103412946918E-06 == Approx( elements[ t21 ].imag() ) );
      CHECK( -5.3265201903090328E-03 == Approx( elements[ t22 ].real() ) );
      CHECK(  3.0309992779761839E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK( -1.9906761702258985E-04 == Approx( elements[ t23 ].real() ) );
      CHECK(  2.5004165149004852E-06 == Approx( elements[ t23 ].imag() ) );
      CHECK( -4.7521348560226043E-04 == Approx( elements[ t31 ].real() ) );
      CHECK(  2.8717398712057162E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK( -1.9906761702258985E-04 == Approx( elements[ t32 ].real() ) );
      CHECK(  2.5004165149004852E-06 == Approx( elements[ t32 ].imag() ) );
      CHECK( -5.4078818278515380E-03 == Approx( elements[ t33 ].real() ) );
      CHECK(  3.0873216863821484E-05 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 1e+3 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -2.3356638158054471E-05 == Approx( elements[ t11 ].real() ) );
      CHECK(  9.7715058458843184E-09 == Approx( elements[ t11 ].imag() ) );
      CHECK( -6.2892983113352601E-05 == Approx( elements[ t12 ].real() ) );
      CHECK(  2.1918820914365044E-08 == Approx( elements[ t12 ].imag() ) );
      CHECK( -6.9651147053680499E-05 == Approx( elements[ t13 ].real() ) );
      CHECK(  3.5260535230245288E-08 == Approx( elements[ t13 ].imag() ) );
      CHECK( -6.2892983113352601E-05 == Approx( elements[ t21 ].real() ) );
      CHECK(  2.1918820914365044E-08 == Approx( elements[ t21 ].imag() ) );
      CHECK( -2.9399211525146107E-04 == Approx( elements[ t22 ].real() ) );
      CHECK(  9.5134378450556152E-08 == Approx( elements[ t22 ].imag() ) );
      CHECK( -1.3971263510252605E-05 == Approx( elements[ t23 ].real() ) );
      CHECK(  1.5068664218306784E-08 == Approx( elements[ t23 ].imag() ) );
      CHECK( -6.9651147053680499E-05 == Approx( elements[ t31 ].real() ) );
      CHECK(  3.5260535230245288E-08 == Approx( elements[ t31 ].imag() ) );
      CHECK( -1.3971263510252605E-05 == Approx( elements[ t32 ].real() ) );
      CHECK(  1.5068664218306784E-08 == Approx( elements[ t32 ].imag() ) );
      CHECK( -4.4949738479358690E-04 == Approx( elements[ t33 ].real() ) );
      CHECK(  2.1641657866935919E-07 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 2e+3 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -1.6226095958907368E-05 == Approx( elements[ t11 ].real() ) );
      CHECK(  3.4220039367265172E-09 == Approx( elements[ t11 ].imag() ) );
      CHECK( -3.6498856824950382E-05 == Approx( elements[ t12 ].real() ) );
      CHECK(  6.3886519761307977E-09 == Approx( elements[ t12 ].imag() ) );
      CHECK( -4.1026702757550493E-05 == Approx( elements[ t13 ].real() ) );
      CHECK(  1.0475988228178095E-08 == Approx( elements[ t13 ].imag() ) );
      CHECK( -3.6498856824950382E-05 == Approx( elements[ t21 ].real() ) );
      CHECK(  6.3886519761307977E-09 == Approx( elements[ t21 ].imag() ) );
      CHECK( -1.4342807016988428E-04 == Approx( elements[ t22 ].real() ) );
      CHECK(  2.3034564822960526E-08 == Approx( elements[ t22 ].imag() ) );
      CHECK( -6.8723140111676477E-06 == Approx( elements[ t23 ].real() ) );
      CHECK(  4.0864233816696672E-09 == Approx( elements[ t23 ].imag() ) );
      CHECK( -4.1026702757550493E-05 == Approx( elements[ t31 ].real() ) );
      CHECK(  1.0475988228178095E-08 == Approx( elements[ t31 ].imag() ) );
      CHECK( -6.8723140111676477E-06 == Approx( elements[ t32 ].real() ) );
      CHECK(  4.0864233816696672E-09 == Approx( elements[ t32 ].imag() ) );
      CHECK( -2.2271662079016285E-04 == Approx( elements[ t33 ].real() ) );
      CHECK(  5.3621475639599757E-08 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 1.541700e+1 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK(  4.6855922882881048E-05 == Approx( elements[ t11 ].real() ) );
      CHECK(  2.5784243657877019E-03 == Approx( elements[ t11 ].imag() ) );
      CHECK(  6.4132524881508777E-04 == Approx( elements[ t12 ].real() ) );
      CHECK(  6.5438089891000262E-05 == Approx( elements[ t12 ].imag() ) );
      CHECK(  5.9795040686584663E-07 == Approx( elements[ t13 ].real() ) );
      CHECK(  4.9399369213564463E-02 == Approx( elements[ t13 ].imag() ) );
      CHECK(  6.4132524881508777E-04 == Approx( elements[ t21 ].real() ) );
      CHECK(  6.5438089891000262E-05 == Approx( elements[ t21 ].imag() ) );
      CHECK(  8.7906104702958170E-03 == Approx( elements[ t22 ].real() ) );
      CHECK(  8.3237160796962133E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK( -2.8267937532975697E-06 == Approx( elements[ t23 ].real() ) );
      CHECK(  1.1394282749259657E-03 == Approx( elements[ t23 ].imag() ) );
      CHECK(  5.9795040686584663E-07 == Approx( elements[ t31 ].real() ) );
      CHECK(  4.9399369213564463E-02 == Approx( elements[ t31 ].imag() ) );
      CHECK( -2.8267937532975697E-06 == Approx( elements[ t32 ].real() ) );
      CHECK(  1.1394282749259657E-03 == Approx( elements[ t32 ].imag() ) );
      CHECK(  9.5700608290543094E-06 == Approx( elements[ t33 ].real() ) );
      CHECK(  0.94658991700278849 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 3.232700e+1 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK(  1.0268471501034247E-04 == Approx( elements[ t11 ].real() ) );
      CHECK(  4.9379169683549618E-03 == Approx( elements[ t11 ].imag() ) );
      CHECK(  1.3355065922058542E-03 == Approx( elements[ t12 ].real() ) );
      CHECK(  1.2147115890292279E-02 == Approx( elements[ t12 ].imag() ) );
      CHECK( -2.2388032699726887E-04 == Approx( elements[ t13 ].real() ) );
      CHECK(  5.9918063277982880E-02 == Approx( elements[ t13 ].imag() ) );
      CHECK(  1.3355065922058542E-03 == Approx( elements[ t21 ].real() ) );
      CHECK(  1.2147115890292279E-02 == Approx( elements[ t21 ].imag() ) );
      CHECK(  1.6647456754967628E-02 == Approx( elements[ t22 ].real() ) );
      CHECK(  3.0081605500153180E-02 == Approx( elements[ t22 ].imag() ) );
      CHECK( -1.8116766863097101E-03 == Approx( elements[ t23 ].real() ) );
      CHECK(  0.14712588608245988 == Approx( elements[ t23 ].imag() ) );
      CHECK( -2.2388032699726887E-04 == Approx( elements[ t31 ].real() ) );
      CHECK(  5.9918063277982880E-02 == Approx( elements[ t31 ].imag() ) );
      CHECK( -1.8116766863097101E-03 == Approx( elements[ t32 ].real() ) );
      CHECK(  0.14712588608245988 == Approx( elements[ t32 ].imag() ) );
      CHECK( -1.1880436411440716E-03 == Approx( elements[ t33 ].real() ) );
      CHECK(  0.72742863159820526 == Approx( elements[ t33 ].imag() ) );

      group2.evaluateTMatrix( 4.753400e+1 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -7.6447370161612620E-05 == Approx( elements[ t11 ].real() ) );
      CHECK(  8.7754624192493379E-03 == Approx( elements[ t11 ].imag() ) );
      CHECK(  4.4048705382639113E-06 == Approx( elements[ t12 ].real() ) );
      CHECK(  9.0882947067441863E-02 == Approx( elements[ t12 ].imag() ) );
      CHECK( -1.1041511969891327E-03 == Approx( elements[ t13 ].real() ) );
      CHECK(  6.2176775301614869E-05 == Approx( elements[ t13 ].imag() ) );
      CHECK(  4.4048705382639113E-06 == Approx( elements[ t21 ].real() ) );
      CHECK(  9.0882947067441863E-02 == Approx( elements[ t21 ].imag() ) );
      CHECK( -5.6430729301813511E-07 == Approx( elements[ t22 ].real() ) );
      CHECK(  0.94136780546751586 == Approx( elements[ t22 ].imag() ) );
      CHECK(  6.5008073784588436E-05 == Approx( elements[ t23 ].real() ) );
      CHECK(  4.4999676874436616E-04 == Approx( elements[ t23 ].imag() ) );
      CHECK( -1.1041511969891327E-03 == Approx( elements[ t31 ].real() ) );
      CHECK(  6.2176775301614869E-05 == Approx( elements[ t31 ].imag() ) );
      CHECK(  6.5008073784588436E-05 == Approx( elements[ t32 ].real() ) );
      CHECK(  4.4999676874436616E-04 == Approx( elements[ t32 ].imag() ) );
      CHECK( -1.5953769305589553E-02 == Approx( elements[ t33 ].real() ) );
      CHECK(  2.6927772660929330E-04 == Approx( elements[ t33 ].imag() ) );
    } // THEN

    THEN( "T matrix elements can be calculated for a single resonance using "
          "the Constant boundary condition" ) {

      std::map< ReactionChannelID, std::complex< double > > elements;
      group3.evaluateTMatrix( 1e-5 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 5.3671960278418717E-08 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.3847864558992168E-09 == Approx( elements[ t11 ].imag() ) );
      CHECK( 4.3622324153697214E-08 == Approx( elements[ t12 ].real() ) );
      CHECK( 1.1254965041247959E-09 == Approx( elements[ t12 ].imag() ) );
      CHECK( 3.6239976480108025E-05 == Approx( elements[ t13 ].real() ) );
      CHECK( 9.3502507326789042E-07 == Approx( elements[ t13 ].imag() ) );
      CHECK( 4.3622324153697214E-08 == Approx( elements[ t21 ].real() ) );
      CHECK( 1.1254965041247959E-09 == Approx( elements[ t21 ].imag() ) );
      CHECK( 3.5454400299505877E-08 == Approx( elements[ t22 ].real() ) );
      CHECK( 9.1475647772318266E-10 == Approx( elements[ t22 ].imag() ) );
      CHECK( 2.9454336922612795E-05 == Approx( elements[ t23 ].real() ) );
      CHECK( 7.5994926636445372E-07 == Approx( elements[ t23 ].imag() ) );
      CHECK( 3.6239976480108025E-05 == Approx( elements[ t31 ].real() ) );
      CHECK( 9.3502507326789042E-07 == Approx( elements[ t31 ].imag() ) );
      CHECK( 2.9454336922612795E-05 == Approx( elements[ t32 ].real() ) );
      CHECK( 7.5994926636445372E-07 == Approx( elements[ t32 ].imag() ) );
      CHECK( 2.4469683769066101E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 6.3134058245241256E-04 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 1e-4 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 1.6972662944694210E-07 == Approx( elements[ t11 ].real() ) );
      CHECK( 4.3791500520677809E-09 == Approx( elements[ t11 ].imag() ) );
      CHECK( 7.7573132643363509E-08 == Approx( elements[ t12 ].real() ) );
      CHECK( 2.0014796084808858E-09 == Approx( elements[ t12 ].imag() ) );
      CHECK( 6.4445179320999683E-05 == Approx( elements[ t13 ].real() ) );
      CHECK( 1.6627627102398558E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK( 7.7573132643363509E-08 == Approx( elements[ t21 ].real() ) );
      CHECK( 2.0014796084808858E-09 == Approx( elements[ t21 ].imag() ) );
      CHECK( 3.5454606785707817E-08 == Approx( elements[ t22 ].real() ) );
      CHECK( 9.1477126281006423E-10 == Approx( elements[ t22 ].imag() ) );
      CHECK( 2.9454508464483770E-05 == Approx( elements[ t23 ].real() ) );
      CHECK( 7.5996154932303602E-07 == Approx( elements[ t23 ].imag() ) );
      CHECK( 6.4445179320999683E-05 == Approx( elements[ t31 ].real() ) );
      CHECK( 1.6627627102398558E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK( 2.9454508464483770E-05 == Approx( elements[ t32 ].real() ) );
      CHECK( 7.5996154932303602E-07 == Approx( elements[ t32 ].imag() ) );
      CHECK( 2.4469826280348800E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 6.3135078672599872E-04 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 1e-3 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 5.3675401113247238E-07 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.3849901146307382E-08 == Approx( elements[ t11 ].imag() ) );
      CHECK( 1.3795474467599038E-07 == Approx( elements[ t12 ].real() ) );
      CHECK( 3.5596558885425563E-09 == Approx( elements[ t12 ].imag() ) );
      CHECK( 1.1460821493055331E-04 == Approx( elements[ t13 ].real() ) );
      CHECK( 2.9572437549072370E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK( 1.3795474467599038E-07 == Approx( elements[ t21 ].real() ) );
      CHECK( 3.5596558885425563E-09 == Approx( elements[ t21 ].imag() ) );
      CHECK( 3.5456673231866484E-08 == Approx( elements[ t22 ].real() ) );
      CHECK( 9.1489100976103619E-10 == Approx( elements[ t22 ].imag() ) );
      CHECK( 2.9456225199243789E-05 == Approx( elements[ t23 ].real() ) );
      CHECK( 7.6006103110835979E-07 == Approx( elements[ t23 ].imag() ) );
      CHECK( 1.1460821493055331E-04 == Approx( elements[ t31 ].real() ) );
      CHECK( 2.9572437549072370E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK( 2.9456225199243789E-05 == Approx( elements[ t32 ].real() ) );
      CHECK( 7.6006103110835979E-07 == Approx( elements[ t32 ].imag() ) );
      CHECK( 2.4471252486506598E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 6.3143343288551165E-04 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 1e-2 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 1.6983553119870929E-06 == Approx( elements[ t11 ].real() ) );
      CHECK( 4.3850353118746114E-08 == Approx( elements[ t11 ].imag() ) );
      CHECK( 2.4546518129987479E-07 == Approx( elements[ t12 ].real() ) );
      CHECK( 6.3377402846067960E-09 == Approx( elements[ t12 ].imag() ) );
      CHECK( 2.0392431099384604E-04 == Approx( elements[ t13 ].real() ) );
      CHECK( 5.2651839008380032E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK( 2.4546518129987479E-07 == Approx( elements[ t21 ].real() ) );
      CHECK( 6.3377402846067960E-09 == Approx( elements[ t21 ].imag() ) );
      CHECK( 3.5477355536447543E-08 == Approx( elements[ t22 ].real() ) );
      CHECK( 9.1600064898807805E-10 == Approx( elements[ t22 ].imag() ) );
      CHECK( 2.9473407370210455E-05 == Approx( elements[ t23 ].real() ) );
      CHECK( 7.6098288248307618E-07 == Approx( elements[ t23 ].imag() ) );
      CHECK( 2.0392431099384604E-04 == Approx( elements[ t31 ].real() ) );
      CHECK( 5.2651839008380032E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK( 2.9473407370210455E-05 == Approx( elements[ t32 ].real() ) );
      CHECK( 7.6098288248307618E-07 == Approx( elements[ t32 ].imag() ) );
      CHECK( 2.4485526862844626E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 6.3219927635639530E-04 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 1e-1 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 5.4021847276018575E-06 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.4032018518734493E-07 == Approx( elements[ t11 ].imag() ) );
      CHECK( 4.3906697669165832E-07 == Approx( elements[ t12 ].real() ) );
      CHECK( 1.1404637676351934E-08 == Approx( elements[ t12 ].imag() ) );
      CHECK( 3.6476224541440953E-04 == Approx( elements[ t13 ].real() ) );
      CHECK( 9.4745937813613308E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK( 4.3906697669165832E-07 == Approx( elements[ t21 ].real() ) );
      CHECK( 1.1404637676351934E-08 == Approx( elements[ t21 ].imag() ) );
      CHECK( 3.5685527197203429E-08 == Approx( elements[ t22 ].real() ) );
      CHECK( 9.2692124340637100E-10 == Approx( elements[ t22 ].imag() ) );
      CHECK( 2.9646349464332672E-05 == Approx( elements[ t23 ].real() ) );
      CHECK( 7.7005534921990856E-07 == Approx( elements[ t23 ].imag() ) );
      CHECK( 3.6476224541440953E-04 == Approx( elements[ t31 ].real() ) );
      CHECK( 9.4745937813613308E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK( 2.9646349464332672E-05 == Approx( elements[ t32 ].real() ) );
      CHECK( 7.7005534921990856E-07 == Approx( elements[ t32 ].imag() ) );
      CHECK( 2.4629201404378163E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 6.3973638006505908E-04 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 1e+0 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 1.8148061780271335E-05 == Approx( elements[ t11 ].real() ) );
      CHECK( 5.0104320258132132E-07 == Approx( elements[ t11 ].imag() ) );
      CHECK( 8.2945263703656714E-07 == Approx( elements[ t12 ].real() ) );
      CHECK( 2.2900054599886333E-08 == Approx( elements[ t12 ].imag() ) );
      CHECK( 6.8908167184441715E-04 == Approx( elements[ t13 ].real() ) );
      CHECK( 1.9024603942904133E-05 == Approx( elements[ t13 ].imag() ) );
      CHECK( 8.2945263703656714E-07 == Approx( elements[ t21 ].real() ) );
      CHECK( 2.2900054599886333E-08 == Approx( elements[ t21 ].imag() ) );
      CHECK( 3.7909925887227664E-08 == Approx( elements[ t22 ].real() ) );
      CHECK( 1.0466412835780584E-09 == Approx( elements[ t22 ].imag() ) );
      CHECK( 3.1494305935538459E-05 == Approx( elements[ t23 ].real() ) );
      CHECK( 8.6951477794573465E-07 == Approx( elements[ t23 ].imag() ) );
      CHECK( 6.8908167184441715E-04 == Approx( elements[ t31 ].real() ) );
      CHECK( 1.9024603942904133E-05 == Approx( elements[ t31 ].imag() ) );
      CHECK( 3.1494305935538459E-05 == Approx( elements[ t32 ].real() ) );
      CHECK( 8.6951477794573465E-07 == Approx( elements[ t32 ].imag() ) );
      CHECK( 2.6164422196759650E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 7.2236396645979800E-04 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 1e+1 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 1.5203097691511610E-04 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.1186916195294557E-05 == Approx( elements[ t11 ].imag() ) );
      CHECK( 3.9074499991002267E-06 == Approx( elements[ t12 ].real() ) );
      CHECK( 2.8752242841696720E-07 == Approx( elements[ t12 ].imag() ) );
      CHECK( 3.2461795379279106E-03 == Approx( elements[ t13 ].real() ) );
      CHECK( 2.3886407351019826E-04 == Approx( elements[ t13 ].imag() ) );
      CHECK( 3.9074499991002267E-06 == Approx( elements[ t21 ].real() ) );
      CHECK( 2.8752242841696720E-07 == Approx( elements[ t21 ].imag() ) );
      CHECK( 1.0042799043508806E-07 == Approx( elements[ t22 ].real() ) );
      CHECK( 7.3898065740013609E-09 == Approx( elements[ t22 ].imag() ) );
      CHECK( 8.3432235258460835E-05 == Approx( elements[ t23 ].real() ) );
      CHECK( 6.1392055932365745E-06 == Approx( elements[ t23 ].imag() ) );
      CHECK( 3.2461795379279106E-03 == Approx( elements[ t31 ].real() ) );
      CHECK( 2.3886407351019826E-04 == Approx( elements[ t31 ].imag() ) );
      CHECK( 8.3432235258460835E-05 == Approx( elements[ t32 ].real() ) );
      CHECK( 6.1392055932365745E-06 == Approx( elements[ t32 ].imag() ) );
      CHECK( 6.9312726960541721E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 5.1002478804555926E-03 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 1e+2 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -3.0955886124614519E-05 == Approx( elements[ t11 ].real() ) );
      CHECK(  1.4653576661553732E-07 == Approx( elements[ t11 ].imag() ) );
      CHECK( -4.4740886351904470E-07 == Approx( elements[ t12 ].real() ) );
      CHECK(  2.1178977252477522E-09 == Approx( elements[ t12 ].imag() ) );
      CHECK( -3.7169240762582824E-04 == Approx( elements[ t13 ].real() ) );
      CHECK(  1.7594790107887367E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK( -4.4740886351904470E-07 == Approx( elements[ t21 ].real() ) );
      CHECK(  2.1178977252477522E-09 == Approx( elements[ t21 ].imag() ) );
      CHECK( -6.4664500427992796E-09 == Approx( elements[ t22 ].real() ) );
      CHECK(  3.0610211269294313E-11 == Approx( elements[ t22 ].imag() ) );
      CHECK( -5.3721116883905769E-06 == Approx( elements[ t23 ].real() ) );
      CHECK(  2.5429945743877634E-08 == Approx( elements[ t23 ].imag() ) );
      CHECK( -3.7169240762582824E-04 == Approx( elements[ t31 ].real() ) );
      CHECK(  1.7594790107887367E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK( -5.3721116883905769E-06 == Approx( elements[ t32 ].real() ) );
      CHECK(  2.5429945743877634E-08 == Approx( elements[ t32 ].imag() ) );
      CHECK( -4.4629717699094058E-03 == Approx( elements[ t33 ].real() ) );
      CHECK(  2.1126353387350176E-05 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 1e+3 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -8.4097609268990072E-06 == Approx( elements[ t11 ].real() ) );
      CHECK(  3.4682701528839583E-09 == Approx( elements[ t11 ].imag() ) );
      CHECK( -6.8351018912904516E-08 == Approx( elements[ t12 ].real() ) );
      CHECK(  2.8188648984846527E-11 == Approx( elements[ t12 ].imag() ) );
      CHECK( -5.6783753865740247E-05 == Approx( elements[ t13 ].real() ) );
      CHECK(  2.3418192314044245E-08 == Approx( elements[ t13 ].imag() ) );
      CHECK( -6.8351018912904516E-08 == Approx( elements[ t21 ].real() ) );
      CHECK(  2.8188648984846527E-11 == Approx( elements[ t21 ].imag() ) );
      CHECK( -5.5552848969690266E-10 == Approx( elements[ t22 ].real() ) );
      CHECK(  2.2910554730869455E-13 == Approx( elements[ t22 ].imag() ) );
      CHECK( -4.6151459811522609E-07 == Approx( elements[ t23 ].real() ) );
      CHECK(  1.9033327102599282E-10 == Approx( elements[ t23 ].imag() ) );
      CHECK( -5.6783753865740247E-05 == Approx( elements[ t31 ].real() ) );
      CHECK(  2.3418192314044245E-08 == Approx( elements[ t31 ].imag() ) );
      CHECK( -4.6151459811522609E-07 == Approx( elements[ t32 ].real() ) );
      CHECK(  1.9033327102599282E-10 == Approx( elements[ t32 ].imag() ) );
      CHECK( -3.8341098291767064E-04 == Approx( elements[ t33 ].real() ) );
      CHECK(  1.5812255305473888E-07 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 2e+3 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -5.9004042760189593E-06 == Approx( elements[ t11 ].real() ) );
      CHECK(  1.2174383646679529E-09 == Approx( elements[ t11 ].imag() ) );
      CHECK( -4.0326043413850607E-08 == Approx( elements[ t12 ].real() ) );
      CHECK(  8.3205268742045395E-12 == Approx( elements[ t12 ].imag() ) );
      CHECK( -3.3501536041607244E-05 == Approx( elements[ t13 ].real() ) );
      CHECK(  6.9124170725259858E-09 == Approx( elements[ t13 ].imag() ) );
      CHECK( -4.0326043413850607E-08 == Approx( elements[ t21 ].real() ) );
      CHECK(  8.3205268742045395E-12 == Approx( elements[ t21 ].imag() ) );
      CHECK( -2.7560650107062913E-10 == Approx( elements[ t22 ].real() ) );
      CHECK(  5.6866260727081882E-14 == Approx( elements[ t22 ].imag() ) );
      CHECK( -2.2896471727121256E-07 == Approx( elements[ t23 ].real() ) );
      CHECK(  4.7242598629089154E-11 == Approx( elements[ t23 ].imag() ) );
      CHECK( -3.3501536041607244E-05 == Approx( elements[ t31 ].real() ) );
      CHECK(  6.9124170725259858E-09 == Approx( elements[ t31 ].imag() ) );
      CHECK( -2.2896471727121256E-07 == Approx( elements[ t32 ].real() ) );
      CHECK(  4.7242598629089154E-11 == Approx( elements[ t32 ].imag() ) );
      CHECK( -1.9021627411340152E-04 == Approx( elements[ t33 ].real() ) );
      CHECK(  3.9247580141423242E-08 == Approx( elements[ t33 ].imag() ) );

      group3.evaluateTMatrix( 1.541700e+1 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -9.9557762280090690E-15 == Approx( elements[ t11 ].real() ) );
      CHECK(  2.5779880786088033E-03 == Approx( elements[ t11 ].imag() ) );
      CHECK( -2.2963416433378427E-16 == Approx( elements[ t12 ].real() ) );
      CHECK(  5.9462378877928035E-05 == Approx( elements[ t12 ].imag() ) );
      CHECK( -1.9077242847409712E-13 == Approx( elements[ t13 ].real() ) );
      CHECK(  4.9399367268805852E-02 == Approx( elements[ t13 ].imag() ) );
      CHECK( -2.2963416433378427E-16 == Approx( elements[ t21 ].real() ) );
      CHECK(  5.9462378877928035E-05 == Approx( elements[ t21 ].imag() ) );
      CHECK( -5.2966085437830095E-18 == Approx( elements[ t22 ].real() ) );
      CHECK(  1.3715247681557664E-06 == Approx( elements[ t22 ].imag() ) );
      CHECK( -4.4002462678224729E-15 == Approx( elements[ t23 ].real() ) );
      CHECK(  1.1394171746724397E-03 == Approx( elements[ t23 ].imag() ) );
      CHECK( -1.9077242847409712E-13 == Approx( elements[ t31 ].real() ) );
      CHECK(  4.9399367268805852E-02 == Approx( elements[ t31 ].imag() ) );
      CHECK( -4.4002462678224729E-15 == Approx( elements[ t32 ].real() ) );
      CHECK(  1.1394171746724397E-03 == Approx( elements[ t32 ].imag() ) );
      CHECK( -3.6555782926814681E-12 == Approx( elements[ t33 ].real() ) );
      CHECK(  0.94658990350152139 == Approx( elements[ t33 ].imag() ) );
    } // THEN

    THEN( "T matrix elements can be calculated for multiple resonances using "
          "the Constant boundary condition" ) {

      std::map< ReactionChannelID, std::complex< double > > elements;
      group4.evaluateTMatrix( 1e-5 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 8.6077679923294300E-08 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.8570343939049694E-09 == Approx( elements[ t11 ].imag() ) );
      CHECK( 1.2888136111385671E-05 == Approx( elements[ t12 ].real() ) );
      CHECK( 9.7894607528028767E-08 == Approx( elements[ t12 ].imag() ) );
      CHECK( 4.0081962349002521E-05 == Approx( elements[ t13 ].real() ) );
      CHECK( 1.1164456176264444E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK( 1.2888136111385671E-05 == Approx( elements[ t21 ].real() ) );
      CHECK( 9.7894607528028767E-08 == Approx( elements[ t21 ].imag() ) );
      CHECK( 5.9174467940555911E-03 == Approx( elements[ t22 ].real() ) );
      CHECK( 3.7061254808378901E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK( 4.3208176704733489E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( 1.4292482833367117E-05 == Approx( elements[ t23 ].imag() ) );
      CHECK( 4.0081962349002521E-05 == Approx( elements[ t31 ].real() ) );
      CHECK( 1.1164456176264444E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK( 4.3208176704733489E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( 1.4292482833367117E-05 == Approx( elements[ t32 ].imag() ) );
      CHECK( 2.6443852465238636E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 7.3345359723414551E-04 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 1e-4 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 2.7220272617802817E-07 == Approx( elements[ t11 ].real() ) );
      CHECK( 5.8725701956782081E-09 == Approx( elements[ t11 ].imag() ) );
      CHECK( 2.2918751872081959E-05 == Approx( elements[ t12 ].real() ) );
      CHECK( 1.7408904822316741E-07 == Approx( elements[ t12 ].imag() ) );
      CHECK( 7.1277322290083349E-05 == Approx( elements[ t13 ].real() ) );
      CHECK( 1.9853875984639007E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK( 2.2918751872081959E-05 == Approx( elements[ t21 ].real() ) );
      CHECK( 1.7408904822316741E-07 == Approx( elements[ t21 ].imag() ) );
      CHECK( 5.9174580639235266E-03 == Approx( elements[ t22 ].real() ) );
      CHECK( 3.7061755644438863E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK( 4.3208301361143982E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( 1.4293712189812671E-05 == Approx( elements[ t23 ].imag() ) );
      CHECK( 7.1277322290083349E-05 == Approx( elements[ t31 ].real() ) );
      CHECK( 1.9853875984639007E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK( 4.3208301361143982E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( 1.4293712189812671E-05 == Approx( elements[ t32 ].imag() ) );
      CHECK( 2.6444000379628831E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 7.3346529717261506E-04 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 1e-3 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 8.6081865742698465E-07 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.8573138778862976E-08 == Approx( elements[ t11 ].imag() ) );
      CHECK( 4.0756742605136367E-05 == Approx( elements[ t12 ].real() ) );
      CHECK( 3.0961750764710787E-07 == Approx( elements[ t12 ].imag() ) );
      CHECK( 1.2675800823402835E-04 == Approx( elements[ t13 ].real() ) );
      CHECK( 3.5310412374359472E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK( 4.0756742605136367E-05 == Approx( elements[ t21 ].real() ) );
      CHECK( 3.0961750764710787E-07 == Approx( elements[ t21 ].imag() ) );
      CHECK( 5.9175708022794109E-03 == Approx( elements[ t22 ].real() ) );
      CHECK( 3.7064308383133275E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK( 4.3209575043568269E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( 1.4298370051136698E-05 == Approx( elements[ t23 ].imag() ) );
      CHECK( 1.2675800823402835E-04 == Approx( elements[ t31 ].real() ) );
      CHECK( 3.5310412374359472E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK( 4.3209575043568269E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( 1.4298370051136698E-05 == Approx( elements[ t32 ].imag() ) );
      CHECK( 2.6445480939713099E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 7.3355857090884812E-04 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 1e-2 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 2.7233522480381398E-06 == Approx( elements[ t11 ].real() ) );
      CHECK( 5.8799731503855438E-08 == Approx( elements[ t11 ].imag() ) );
      CHECK( 7.2491079343976574E-05 == Approx( elements[ t12 ].real() ) );
      CHECK( 5.5098032781778119E-07 == Approx( elements[ t12 ].imag() ) );
      CHECK( 2.2553600034014542E-04 == Approx( elements[ t13 ].real() ) );
      CHECK( 6.2865891242128442E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK( 7.2491079343976574E-05 == Approx( elements[ t21 ].real() ) );
      CHECK( 5.5098032781778119E-07 == Approx( elements[ t21 ].imag() ) );
      CHECK( 5.9186985408499057E-03 == Approx( elements[ t22 ].real() ) );
      CHECK( 3.7082076123237578E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK( 4.3222401908481596E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( 1.4320818280330810E-05 == Approx( elements[ t23 ].imag() ) );
      CHECK( 2.2553600034014542E-04 == Approx( elements[ t31 ].real() ) );
      CHECK( 6.2865891242128442E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK( 4.3222401908481596E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( 1.4320818280330810E-05 == Approx( elements[ t32 ].imag() ) );
      CHECK( 2.6460300033077304E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 7.3441706653197981E-04 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 1e-1 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 8.6503127273711784E-06 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.8794672394205753E-07 == Approx( elements[ t11 ].imag() ) );
      CHECK( 1.2916257209232690E-04 == Approx( elements[ t12 ].real() ) );
      CHECK( 9.8518792030730516E-07 == Approx( elements[ t12 ].imag() ) );
      CHECK( 4.0330031176166552E-04 == Approx( elements[ t13 ].real() ) );
      CHECK( 1.1307201349134586E-05 == Approx( elements[ t13 ].imag() ) );
      CHECK( 1.2916257209232690E-04 == Approx( elements[ t21 ].real() ) );
      CHECK( 9.8518792030730516E-07 == Approx( elements[ t21 ].imag() ) );
      CHECK( 5.9300000563923538E-03 == Approx( elements[ t22 ].real() ) );
      CHECK( 3.7235694751397835E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK( 4.3351388772216185E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( 1.4470029756190969E-05 == Approx( elements[ t23 ].imag() ) );
      CHECK( 4.0330031176166552E-04 == Approx( elements[ t31 ].real() ) );
      CHECK( 1.1307201349134586E-05 == Approx( elements[ t31 ].imag() ) );
      CHECK( 4.3351388772216185E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( 1.4470029756190969E-05 == Approx( elements[ t32 ].imag() ) );
      CHECK( 2.6609441240777808E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 7.4284391249763656E-04 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 1e+0 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 2.8639910507909060E-05 == Approx( elements[ t11 ].real() ) );
      CHECK( 6.6284973829173364E-07 == Approx( elements[ t11 ].imag() ) );
      CHECK( 2.3429141040536179E-04 == Approx( elements[ t12 ].real() ) );
      CHECK( 1.8428769093161273E-06 == Approx( elements[ t12 ].imag() ) );
      CHECK( 7.5956489245821639E-04 == Approx( elements[ t13 ].real() ) );
      CHECK( 2.2581707931633565E-05 == Approx( elements[ t13 ].imag() ) );
      CHECK( 2.3429141040536179E-04 == Approx( elements[ t21 ].real() ) );
      CHECK( 1.8428769093161273E-06 == Approx( elements[ t21 ].imag() ) );
      CHECK( 6.0454438704396537E-03 == Approx( elements[ t22 ].real() ) );
      CHECK( 3.8742163879110548E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK( 4.4688779003336233E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( 1.5816695069991596E-05 == Approx( elements[ t23 ].imag() ) );
      CHECK( 7.5956489245821639E-04 == Approx( elements[ t31 ].real() ) );
      CHECK( 2.2581707931633565E-05 == Approx( elements[ t31 ].imag() ) );
      CHECK( 4.4688779003336233E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( 1.5816695069991596E-05 == Approx( elements[ t32 ].imag() ) );
      CHECK( 2.8201016434418785E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 8.3491814606862077E-04 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 1e+1 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( 1.9431416180186621E-04 == Approx( elements[ t11 ].real() ) );
      CHECK( 1.2666815936613832E-05 == Approx( elements[ t11 ].imag() ) );
      CHECK( 5.2259368632340708E-04 == Approx( elements[ t12 ].real() ) );
      CHECK( 6.5435669882954769E-06 == Approx( elements[ t12 ].imag() ) );
      CHECK( 3.4199812581358797E-03 == Approx( elements[ t13 ].real() ) );
      CHECK( 2.6145017724662367E-04 == Approx( elements[ t13 ].imag() ) );
      CHECK( 5.2259368632340708E-04 == Approx( elements[ t21 ].real() ) );
      CHECK( 6.5435669882954769E-06 == Approx( elements[ t21 ].imag() ) );
      CHECK( 7.5085348211873038E-03 == Approx( elements[ t22 ].real() ) );
      CHECK( 6.0100811547122537E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK( 6.6260767326451852E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( 5.5721401343776834E-05 == Approx( elements[ t23 ].imag() ) );
      CHECK( 3.4199812581358797E-03 == Approx( elements[ t31 ].real() ) );
      CHECK( 2.6145017724662367E-04 == Approx( elements[ t31 ].imag() ) );
      CHECK( 6.6260767326451852E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( 5.5721401343776834E-05 == Approx( elements[ t32 ].imag() ) );
      CHECK( 7.2131730438568478E-02 == Approx( elements[ t33 ].real() ) );
      CHECK( 5.5074649724989677E-03 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 1e+2 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -1.1371918230606257E-04 == Approx( elements[ t11 ].real() ) );
      CHECK(  6.7440028937519096E-07 == Approx( elements[ t11 ].imag() ) );
      CHECK( -6.3615465189479985E-04 == Approx( elements[ t12 ].real() ) );
      CHECK(  3.7343103412946918E-06 == Approx( elements[ t12 ].imag() ) );
      CHECK( -4.7521348560226043E-04 == Approx( elements[ t13 ].real() ) );
      CHECK(  2.8717398712057162E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK( -6.3615465189479985E-04 == Approx( elements[ t21 ].real() ) );
      CHECK(  3.7343103412946918E-06 == Approx( elements[ t21 ].imag() ) );
      CHECK( -5.3265201903090328E-03 == Approx( elements[ t22 ].real() ) );
      CHECK(  3.0309992779761839E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK( -1.9906761702258985E-04 == Approx( elements[ t23 ].real() ) );
      CHECK(  2.5004165149004852E-06 == Approx( elements[ t23 ].imag() ) );
      CHECK( -4.7521348560226043E-04 == Approx( elements[ t31 ].real() ) );
      CHECK(  2.8717398712057162E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK( -1.9906761702258985E-04 == Approx( elements[ t32 ].real() ) );
      CHECK(  2.5004165149004852E-06 == Approx( elements[ t32 ].imag() ) );
      CHECK( -5.4078818278515380E-03 == Approx( elements[ t33 ].real() ) );
      CHECK(  3.0873216863821484E-05 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 1e+3 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -2.3356638158054471E-05 == Approx( elements[ t11 ].real() ) );
      CHECK(  9.7715058458843184E-09 == Approx( elements[ t11 ].imag() ) );
      CHECK( -6.2892983113352601E-05 == Approx( elements[ t12 ].real() ) );
      CHECK(  2.1918820914365044E-08 == Approx( elements[ t12 ].imag() ) );
      CHECK( -6.9651147053680499E-05 == Approx( elements[ t13 ].real() ) );
      CHECK(  3.5260535230245288E-08 == Approx( elements[ t13 ].imag() ) );
      CHECK( -6.2892983113352601E-05 == Approx( elements[ t21 ].real() ) );
      CHECK(  2.1918820914365044E-08 == Approx( elements[ t21 ].imag() ) );
      CHECK( -2.9399211525146107E-04 == Approx( elements[ t22 ].real() ) );
      CHECK(  9.5134378450556152E-08 == Approx( elements[ t22 ].imag() ) );
      CHECK( -1.3971263510252605E-05 == Approx( elements[ t23 ].real() ) );
      CHECK(  1.5068664218306784E-08 == Approx( elements[ t23 ].imag() ) );
      CHECK( -6.9651147053680499E-05 == Approx( elements[ t31 ].real() ) );
      CHECK(  3.5260535230245288E-08 == Approx( elements[ t31 ].imag() ) );
      CHECK( -1.3971263510252605E-05 == Approx( elements[ t32 ].real() ) );
      CHECK(  1.5068664218306784E-08 == Approx( elements[ t32 ].imag() ) );
      CHECK( -4.4949738479358690E-04 == Approx( elements[ t33 ].real() ) );
      CHECK(  2.1641657866935919E-07 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 2e+3 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -1.6226095958907368E-05 == Approx( elements[ t11 ].real() ) );
      CHECK(  3.4220039367265172E-09 == Approx( elements[ t11 ].imag() ) );
      CHECK( -3.6498856824950382E-05 == Approx( elements[ t12 ].real() ) );
      CHECK(  6.3886519761307977E-09 == Approx( elements[ t12 ].imag() ) );
      CHECK( -4.1026702757550493E-05 == Approx( elements[ t13 ].real() ) );
      CHECK(  1.0475988228178095E-08 == Approx( elements[ t13 ].imag() ) );
      CHECK( -3.6498856824950382E-05 == Approx( elements[ t21 ].real() ) );
      CHECK(  6.3886519761307977E-09 == Approx( elements[ t21 ].imag() ) );
      CHECK( -1.4342807016988428E-04 == Approx( elements[ t22 ].real() ) );
      CHECK(  2.3034564822960526E-08 == Approx( elements[ t22 ].imag() ) );
      CHECK( -6.8723140111676477E-06 == Approx( elements[ t23 ].real() ) );
      CHECK(  4.0864233816696672E-09 == Approx( elements[ t23 ].imag() ) );
      CHECK( -4.1026702757550493E-05 == Approx( elements[ t31 ].real() ) );
      CHECK(  1.0475988228178095E-08 == Approx( elements[ t31 ].imag() ) );
      CHECK( -6.8723140111676477E-06 == Approx( elements[ t32 ].real() ) );
      CHECK(  4.0864233816696672E-09 == Approx( elements[ t32 ].imag() ) );
      CHECK( -2.2271662079016285E-04 == Approx( elements[ t33 ].real() ) );
      CHECK(  5.3621475639599757E-08 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 1.541700e+1 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK(  4.6855922882881048E-05 == Approx( elements[ t11 ].real() ) );
      CHECK(  2.5784243657877019E-03 == Approx( elements[ t11 ].imag() ) );
      CHECK(  6.4132524881508777E-04 == Approx( elements[ t12 ].real() ) );
      CHECK(  6.5438089891000262E-05 == Approx( elements[ t12 ].imag() ) );
      CHECK(  5.9795040686584663E-07 == Approx( elements[ t13 ].real() ) );
      CHECK(  4.9399369213564463E-02 == Approx( elements[ t13 ].imag() ) );
      CHECK(  6.4132524881508777E-04 == Approx( elements[ t21 ].real() ) );
      CHECK(  6.5438089891000262E-05 == Approx( elements[ t21 ].imag() ) );
      CHECK(  8.7906104702958170E-03 == Approx( elements[ t22 ].real() ) );
      CHECK(  8.3237160796962133E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK( -2.8267937532975697E-06 == Approx( elements[ t23 ].real() ) );
      CHECK(  1.1394282749259657E-03 == Approx( elements[ t23 ].imag() ) );
      CHECK(  5.9795040686584663E-07 == Approx( elements[ t31 ].real() ) );
      CHECK(  4.9399369213564463E-02 == Approx( elements[ t31 ].imag() ) );
      CHECK( -2.8267937532975697E-06 == Approx( elements[ t32 ].real() ) );
      CHECK(  1.1394282749259657E-03 == Approx( elements[ t32 ].imag() ) );
      CHECK(  9.5700608290543094E-06 == Approx( elements[ t33 ].real() ) );
      CHECK(  0.94658991700278849 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 3.232700e+1 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK(  1.0268471501034247E-04 == Approx( elements[ t11 ].real() ) );
      CHECK(  4.9379169683549618E-03 == Approx( elements[ t11 ].imag() ) );
      CHECK(  1.3355065922058542E-03 == Approx( elements[ t12 ].real() ) );
      CHECK(  1.2147115890292279E-02 == Approx( elements[ t12 ].imag() ) );
      CHECK( -2.2388032699726887E-04 == Approx( elements[ t13 ].real() ) );
      CHECK(  5.9918063277982880E-02 == Approx( elements[ t13 ].imag() ) );
      CHECK(  1.3355065922058542E-03 == Approx( elements[ t21 ].real() ) );
      CHECK(  1.2147115890292279E-02 == Approx( elements[ t21 ].imag() ) );
      CHECK(  1.6647456754967628E-02 == Approx( elements[ t22 ].real() ) );
      CHECK(  3.0081605500153180E-02 == Approx( elements[ t22 ].imag() ) );
      CHECK( -1.8116766863097101E-03 == Approx( elements[ t23 ].real() ) );
      CHECK(  0.14712588608245988 == Approx( elements[ t23 ].imag() ) );
      CHECK( -2.2388032699726887E-04 == Approx( elements[ t31 ].real() ) );
      CHECK(  5.9918063277982880E-02 == Approx( elements[ t31 ].imag() ) );
      CHECK( -1.8116766863097101E-03 == Approx( elements[ t32 ].real() ) );
      CHECK(  0.14712588608245988 == Approx( elements[ t32 ].imag() ) );
      CHECK( -1.1880436411440716E-03 == Approx( elements[ t33 ].real() ) );
      CHECK(  0.72742863159820526 == Approx( elements[ t33 ].imag() ) );

      group4.evaluateTMatrix( 4.753400e+1 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -7.6447370161612620E-05 == Approx( elements[ t11 ].real() ) );
      CHECK(  8.7754624192493379E-03 == Approx( elements[ t11 ].imag() ) );
      CHECK(  4.4048705382639113E-06 == Approx( elements[ t12 ].real() ) );
      CHECK(  9.0882947067441863E-02 == Approx( elements[ t12 ].imag() ) );
      CHECK( -1.1041511969891327E-03 == Approx( elements[ t13 ].real() ) );
      CHECK(  6.2176775301614869E-05 == Approx( elements[ t13 ].imag() ) );
      CHECK(  4.4048705382639113E-06 == Approx( elements[ t21 ].real() ) );
      CHECK(  9.0882947067441863E-02 == Approx( elements[ t21 ].imag() ) );
      CHECK( -5.6430729301813511E-07 == Approx( elements[ t22 ].real() ) );
      CHECK(  0.94136780546751586 == Approx( elements[ t22 ].imag() ) );
      CHECK(  6.5008073784588436E-05 == Approx( elements[ t23 ].real() ) );
      CHECK(  4.4999676874436616E-04 == Approx( elements[ t23 ].imag() ) );
      CHECK( -1.1041511969891327E-03 == Approx( elements[ t31 ].real() ) );
      CHECK(  6.2176775301614869E-05 == Approx( elements[ t31 ].imag() ) );
      CHECK(  6.5008073784588436E-05 == Approx( elements[ t32 ].real() ) );
      CHECK(  4.4999676874436616E-04 == Approx( elements[ t32 ].imag() ) );
      CHECK( -1.5953769305589553E-02 == Approx( elements[ t33 ].real() ) );
      CHECK(  2.6927772660929330E-04 == Approx( elements[ t33 ].imag() ) );
    } // THEN
  } // GIVEN

  GIVEN( "valid data for a SpinGroup with a resonance at a negative energy" ) {

    // test based on Pu239 ENDF/B-VIII.0 LRF3 resonance evaluation (some of the
    // widths used in this test were negative in the original evaluation)
    // data given in Gamma = 2 gamma^2 P(Er) so conversion is required
    // T matrix values extracted from the samm.f90 routines in NJOY2016.43
    // (note: LRF7 in NJOY2016 doesn't add potential scattering for missing J
    // values)

    // particles
    Particle photon( ParticleID( "g" ), 0.0 * daltons,
                     0.0 * elementary, 1., +1);
    Particle neutron( ParticleID( "n" ), neutronMass,
                      0.0 * elementary, 0.5, +1);
    Particle pu240( ParticleID( "Pu240_e0" ), 2.379916e+2 * neutronMass,
                    94.0 * elementary, 0.5, +1);
    Particle pu239( ParticleID( "Pu239_e0" ), 2.369986e+2 * neutronMass,
                    94.0 * elementary, 0.5, +1);

    // particle pairs
    ParticlePair in( neutron, pu239 );
    ParticlePair out1( photon, pu240 );
    ParticlePair out2( neutron, pu239, ParticlePairID( "fission" ) );

    // channels
    Channel< Photon > capture( in, out1, 0. * electronVolt, { 0, 0.0, 0.0, +1 },
                               { 0.0 * rootBarn },
                               0.0 );
    Channel< Neutron > elastic( in, in, 0. * electronVolt, { 0, 0.5, 0.0, +1 },
                                { 9.410000e-1 * rootBarn },
                                0.0 );
    Channel< Fission > fission1( in, out2, "fission1", 0. * electronVolt,
                                 { 0, 0.0, 0.0, +1 },
                                 { 9.410000e-1 * rootBarn },
                                 0.0 );
    Channel< Fission > fission2( in, out2, "fission2", 0. * electronVolt,
                                 { 0, 0.0, 0.0, +1 },
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
        group( { elastic, fission1, fission2 }, std::move( table ) );

    ReactionChannelID t11( "n,Pu239{0,1/2,0+}->n,Pu239{0,1/2,0+}" );
    ReactionChannelID t12( "n,Pu239{0,1/2,0+}->fission1{0,0,0+}" );
    ReactionChannelID t13( "n,Pu239{0,1/2,0+}->fission2{0,0,0+}" );
    ReactionChannelID t21( "fission1{0,0,0+}->n,Pu239{0,1/2,0+}" );
    ReactionChannelID t22( "fission1{0,0,0+}->fission1{0,0,0+}" );
    ReactionChannelID t23( "fission1{0,0,0+}->fission2{0,0,0+}" );
    ReactionChannelID t31( "fission2{0,0,0+}->n,Pu239{0,1/2,0+}" );
    ReactionChannelID t32( "fission2{0,0,0+}->fission1{0,0,0+}" );
    ReactionChannelID t33( "fission2{0,0,0+}->fission2{0,0,0+}" );

    THEN( "T matrix elements can be calculated" ) {

      std::map< ReactionChannelID, std::complex< double > > elements;
      group.evaluateTMatrix( 1e-5 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -2.1263213571644141E-08 == Approx( elements[ t11 ].real() ) );
      CHECK(  1.2971372284823744E-09 == Approx( elements[ t11 ].imag() ) );
      CHECK(  1.2801222616473466E-05 == Approx( elements[ t12 ].real() ) );
      CHECK(  6.7941497903637647E-08 == Approx( elements[ t12 ].imag() ) );
      CHECK( -3.2396429913820167E-05 == Approx( elements[ t13 ].real() ) );
      CHECK(  7.8409475874204806E-07 == Approx( elements[ t13 ].imag() ) );
      CHECK(  1.2801222616473466E-05 == Approx( elements[ t21 ].real() ) );
      CHECK(  6.7941497903637647E-08 == Approx( elements[ t21 ].imag() ) );
      CHECK(  5.9173845275415629E-03 == Approx( elements[ t22 ].real() ) );
      CHECK(  3.7013073766943528E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK(  3.7338374143853825E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( -5.8952953266736964E-06 == Approx( elements[ t23 ].imag() ) );
      CHECK( -3.2396429913820167E-05 == Approx( elements[ t31 ].real() ) );
      CHECK(  7.8409475874204806E-07 == Approx( elements[ t31 ].imag() ) );
      CHECK(  3.7338374143853825E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( -5.8952953266736964E-06 == Approx( elements[ t32 ].imag() ) );
      CHECK( -2.2494768086750154E-02 == Approx( elements[ t33 ].real() ) );
      CHECK(  5.3990278696026061E-04 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e-4 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -6.7238980019747505E-08 == Approx( elements[ t11 ].real() ) );
      CHECK(  4.1018652337798205E-09 == Approx( elements[ t11 ].imag() ) );
      CHECK(  2.2764195478870796E-05 == Approx( elements[ t12 ].real() ) );
      CHECK(  1.2081854447813631E-07 == Approx( elements[ t12 ].imag() ) );
      CHECK( -5.7609509202257288E-05 == Approx( elements[ t13 ].real() ) );
      CHECK(  1.3943238526046031E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK(  2.2764195478870796E-05 == Approx( elements[ t21 ].real() ) );
      CHECK(  1.2081854447813631E-07 == Approx( elements[ t21 ].imag() ) );
      CHECK(  5.9173957991666181E-03 == Approx( elements[ t22 ].real() ) );
      CHECK(  3.7013569634782119E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK(  3.7338501782878684E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( -5.8961530002853961E-06 == Approx( elements[ t23 ].imag() ) );
      CHECK( -5.7609509202257288E-05 == Approx( elements[ t31 ].real() ) );
      CHECK(  1.3943238526046031E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK(  3.7338501782878684E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( -5.8961530002853961E-06 == Approx( elements[ t32 ].imag() ) );
      CHECK( -2.2494619781134262E-02 == Approx( elements[ t33 ].real() ) );
      CHECK(  5.3989801277359297E-04 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e-3 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -2.1259022133504974E-07 == Approx( elements[ t11 ].real() ) );
      CHECK(  1.2969815423711925E-08 == Approx( elements[ t11 ].imag() ) );
      CHECK(  4.0481898864244864E-05 == Approx( elements[ t12 ].real() ) );
      CHECK(  2.1485435148825168E-07 == Approx( elements[ t12 ].imag() ) );
      CHECK( -1.0243878171068415E-04 == Approx( elements[ t13 ].real() ) );
      CHECK(  2.4791858563276500E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK(  4.0481898864244864E-05 == Approx( elements[ t21 ].real() ) );
      CHECK(  2.1485435148825168E-07 == Approx( elements[ t21 ].imag() ) );
      CHECK(  5.9175085435045054E-03 == Approx( elements[ t22 ].real() ) );
      CHECK(  3.7016105628664716E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK(  3.7339789782970438E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( -5.8985972415680078E-06 == Approx( elements[ t23 ].imag() ) );
      CHECK( -1.0243878171068415E-04 == Approx( elements[ t31 ].real() ) );
      CHECK(  2.4791858563276500E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK(  3.7339789782970438E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( -5.8985972415680078E-06 == Approx( elements[ t32 ].imag() ) );
      CHECK( -2.2493137565948215E-02 == Approx( elements[ t33 ].real() ) );
      CHECK(  5.3983476840927208E-04 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e-2 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -6.7106504705850411E-07 == Approx( elements[ t11 ].real() ) );
      CHECK(  4.0968590226101187E-08 == Approx( elements[ t11 ].imag() ) );
      CHECK(  7.2002338323246527E-05 == Approx( elements[ t12 ].real() ) );
      CHECK(  3.8223560913550099E-07 == Approx( elements[ t12 ].imag() ) );
      CHECK( -1.8203998944443065E-04 == Approx( elements[ t13 ].real() ) );
      CHECK(  4.4029696690378768E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK(  7.2002338323246527E-05 == Approx( elements[ t21 ].real() ) );
      CHECK(  3.8223560913550099E-07 == Approx( elements[ t21 ].imag() ) );
      CHECK(  5.9186363052359282E-03 == Approx( elements[ t22 ].real() ) );
      CHECK(  3.7033810048139998E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK(  3.7352708751740404E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( -5.9036440015630715E-06 == Approx( elements[ t23 ].imag() ) );
      CHECK( -1.8203998944443065E-04 == Approx( elements[ t31 ].real() ) );
      CHECK(  4.4029696690378768E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK(  3.7352708751740404E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( -5.9036440015630715E-06 == Approx( elements[ t32 ].imag() ) );
      CHECK( -2.2478326761901046E-02 == Approx( elements[ t33 ].real() ) );
      CHECK(  5.3915393248015278E-04 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e-1 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -2.0841956999526724E-06 == Approx( elements[ t11 ].real() ) );
      CHECK(  1.2812120650339705E-07 == Approx( elements[ t11 ].imag() ) );
      CHECK(  1.2829353810576006E-04 == Approx( elements[ t12 ].real() ) );
      CHECK(  6.8306086753656647E-07 == Approx( elements[ t12 ].imag() ) );
      CHECK( -3.2151211300204835E-04 == Approx( elements[ t13 ].real() ) );
      CHECK(  7.7280408405426308E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK(  1.2829353810576006E-04 == Approx( elements[ t21 ].real() ) );
      CHECK(  6.8306086753656647E-07 == Approx( elements[ t21 ].imag() ) );
      CHECK(  5.9299379347602758E-03 == Approx( elements[ t22 ].real() ) );
      CHECK(  3.7187122345964999E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK(  3.7482247312349001E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( -5.8927666460044360E-06 == Approx( elements[ t23 ].imag() ) );
      CHECK( -3.2151211300204835E-04 == Approx( elements[ t31 ].real() ) );
      CHECK(  7.7280408405426308E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK(  3.7482247312349001E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( -5.8927666460044360E-06 == Approx( elements[ t32 ].imag() ) );
      CHECK( -2.2331117399884647E-02 == Approx( elements[ t33 ].real() ) );
      CHECK(  5.3225711342881816E-04 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e+0 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -5.4462576894535023E-06 == Approx( elements[ t11 ].real() ) );
      CHECK(  3.6410337482659483E-07 == Approx( elements[ t11 ].imag() ) );
      CHECK(  2.3274167939792254E-04 == Approx( elements[ t12 ].real() ) );
      CHECK(  1.2778852474274546E-06 == Approx( elements[ t12 ].imag() ) );
      CHECK( -5.3469862330085462E-04 == Approx( elements[ t13 ].real() ) );
      CHECK(  1.2091349386787632E-05 == Approx( elements[ t13 ].imag() ) );
      CHECK(  2.3274167939792254E-04 == Approx( elements[ t21 ].real() ) );
      CHECK(  1.2778852474274546E-06 == Approx( elements[ t21 ].imag() ) );
      CHECK(  6.0453823311809726E-03 == Approx( elements[ t22 ].real() ) );
      CHECK(  3.8691328433044154E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK(  3.8803017114425169E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( -5.5973232137600748E-06 == Approx( elements[ t23 ].imag() ) );
      CHECK( -5.3469862330085462E-04 == Approx( elements[ t31 ].real() ) );
      CHECK(  1.2091349386787632E-05 == Approx( elements[ t31 ].imag() ) );
      CHECK(  3.8803017114425169E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( -5.5973232137600748E-06 == Approx( elements[ t32 ].imag() ) );
      CHECK( -2.0942592513586000E-02 == Approx( elements[ t33 ].real() ) );
      CHECK(  4.6898959315351633E-04 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e+1 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK(  9.8349840297660857E-06 == Approx( elements[ t11 ].real() ) );
      CHECK(  5.8666424348490828E-07 == Approx( elements[ t11 ].imag() ) );
      CHECK(  5.1802814213470138E-04 == Approx( elements[ t12 ].real() ) );
      CHECK(  3.8246542457674067E-06 == Approx( elements[ t12 ].imag() ) );
      CHECK( -5.1930424677988933E-04 == Approx( elements[ t13 ].real() ) );
      CHECK(  7.2326367453940488E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK(  5.1802814213470138E-04 == Approx( elements[ t21 ].real() ) );
      CHECK(  3.8246542457674067E-06 == Approx( elements[ t21 ].imag() ) );
      CHECK(  7.5084534405060466E-03 == Approx( elements[ t22 ].real() ) );
      CHECK(  5.9971563672885391E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK(  5.6506828686182461E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( -2.2413699337732195E-06 == Approx( elements[ t23 ].imag() ) );
      CHECK( -5.1930424677988933E-04 == Approx( elements[ t31 ].real() ) );
      CHECK(  7.2326367453940488E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK(  5.6506828686182461E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( -2.2413699337732195E-06 == Approx( elements[ t32 ].imag() ) );
      CHECK( -1.1985927576449064E-02 == Approx( elements[ t33 ].real() ) );
      CHECK(  1.5880238733464967E-04 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e+2 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -1.0544977684396769E-04 == Approx( elements[ t11 ].real() ) );
      CHECK(  5.8448538114818004E-07 == Approx( elements[ t11 ].imag() ) );
      CHECK( -6.3603553777656727E-04 == Approx( elements[ t12 ].real() ) );
      CHECK(  3.7080441601114399E-06 == Approx( elements[ t12 ].imag() ) );
      CHECK( -3.7592077890816607E-04 == Approx( elements[ t13 ].real() ) );
      CHECK(  1.8299654650210189E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK( -6.3603553777656727E-04 == Approx( elements[ t21 ].real() ) );
      CHECK(  3.7080441601114399E-06 == Approx( elements[ t21 ].imag() ) );
      CHECK( -5.3265185499620166E-03 == Approx( elements[ t22 ].real() ) );
      CHECK(  3.0309253929260578E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK( -1.9763727134302524E-04 == Approx( elements[ t23 ].real() ) );
      CHECK(  2.1855787825610621E-06 == Approx( elements[ t23 ].imag() ) );
      CHECK( -3.7592077890816607E-04 == Approx( elements[ t31 ].real() ) );
      CHECK(  1.8299654650210189E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK( -1.9763727134302524E-04 == Approx( elements[ t32 ].real() ) );
      CHECK(  2.1855787825610621E-06 == Approx( elements[ t32 ].imag() ) );
      CHECK( -4.2156509939678209E-03 == Approx( elements[ t33 ].real() ) );
      CHECK(  1.8818919181566614E-05 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e+3 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -2.3101268773318971E-05 == Approx( elements[ t11 ].real() ) );
      CHECK(  9.5118011008826436E-09 == Approx( elements[ t11 ].imag() ) );
      CHECK( -6.2890907629854781E-05 == Approx( elements[ t12 ].real() ) );
      CHECK(  2.1876974201514174E-08 == Approx( elements[ t12 ].imag() ) );
      CHECK( -6.7926861065473792E-05 == Approx( elements[ t13 ].real() ) );
      CHECK(  3.3566173559980896E-08 == Approx( elements[ t13 ].imag() ) );
      CHECK( -6.2890907629854781E-05 == Approx( elements[ t21 ].real() ) );
      CHECK(  2.1876974201514174E-08 == Approx( elements[ t21 ].imag() ) );
      CHECK( -2.9399209838940532E-04 == Approx( elements[ t22 ].real() ) );
      CHECK(  9.5133715390776343E-08 == Approx( elements[ t22 ].imag() ) );
      CHECK( -1.3957249577180666E-05 == Approx( elements[ t23 ].real() ) );
      CHECK(  1.4786591112844456E-08 == Approx( elements[ t23 ].imag() ) );
      CHECK( -6.7926861065473792E-05 == Approx( elements[ t31 ].real() ) );
      CHECK(  3.3566173559980896E-08 == Approx( elements[ t31 ].imag() ) );
      CHECK( -1.3957249577180666E-05 == Approx( elements[ t32 ].real() ) );
      CHECK(  1.4786591112844456E-08 == Approx( elements[ t32 ].imag() ) );
      CHECK( -4.3785479040545707E-04 == Approx( elements[ t33 ].real() ) );
      CHECK(  2.0537574395780616E-07 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 2e+3 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -1.6135825291599936E-05 == Approx( elements[ t11 ].real() ) );
      CHECK(  3.3754144197369945E-09 == Approx( elements[ t11 ].imag() ) );
      CHECK( -3.6498239878646625E-05 == Approx( elements[ t12 ].real() ) );
      CHECK(  6.3815818260073107E-09 == Approx( elements[ t12 ].imag() ) );
      CHECK( -4.0514160580584947E-05 == Approx( elements[ t13 ].real() ) );
      CHECK(  1.0221438614421095E-08 == Approx( elements[ t13 ].imag() ) );
      CHECK( -3.6498239878646625E-05 == Approx( elements[ t21 ].real() ) );
      CHECK(  6.3815818260073107E-09 == Approx( elements[ t21 ].imag() ) );
      CHECK( -1.4342806595392830E-04 == Approx( elements[ t22 ].real() ) );
      CHECK(  2.3034470358314618E-08 == Approx( elements[ t22 ].imag() ) );
      CHECK( -6.8688110895830874E-06 == Approx( elements[ t23 ].real() ) );
      CHECK(  4.0463484131077039E-09 == Approx( elements[ t23 ].imag() ) );
      CHECK( -4.0514160580584947E-05 == Approx( elements[ t31 ].real() ) );
      CHECK(  1.0221438614421095E-08 == Approx( elements[ t31 ].imag() ) );
      CHECK( -6.8688110895830874E-06 == Approx( elements[ t32 ].real() ) );
      CHECK(  4.0463484131077039E-09 == Approx( elements[ t32 ].imag() ) );
      CHECK( -2.1980648961572108E-04 == Approx( elements[ t33 ].real() ) );
      CHECK(  5.2232838285882173E-08 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 3.232700e+1 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK(  1.0277826749300980E-04 == Approx( elements[ t11 ].real() ) );
      CHECK(  4.9379012722419050E-03 == Approx( elements[ t11 ].imag() ) );
      CHECK(  1.3300936419776697E-03 == Approx( elements[ t12 ].real() ) );
      CHECK(  1.2147497328949593E-02 == Approx( elements[ t12 ].imag() ) );
      CHECK( -2.1398288784814625E-04 == Approx( elements[ t13 ].real() ) );
      CHECK(  5.9917205271759158E-02 == Approx( elements[ t13 ].imag() ) );
      CHECK(  1.3300936419776697E-03 == Approx( elements[ t21 ].real() ) );
      CHECK(  1.2147497328949593E-02 == Approx( elements[ t21 ].imag() ) );
      CHECK(  1.6957765703829918E-02 == Approx( elements[ t22 ].real() ) );
      CHECK(  3.0089528618649636E-02 == Approx( elements[ t22 ].imag() ) );
      CHECK( -2.3799477972216740E-03 == Approx( elements[ t23 ].real() ) );
      CHECK(  0.14712054109643591 == Approx( elements[ t23 ].imag() ) );
      CHECK( -2.1398288784814625E-04 == Approx( elements[ t31 ].real() ) );
      CHECK(  5.9917205271759158E-02 == Approx( elements[ t31 ].imag() ) );
      CHECK( -2.3799477972216740E-03 == Approx( elements[ t32 ].real() ) );
      CHECK(  0.14712054109643591 == Approx( elements[ t32 ].imag() ) );
      CHECK( -1.4763497734517853E-04 == Approx( elements[ t33 ].real() ) );
      CHECK(  0.72742164346961768 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 4.753400e+1 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -4.9533155268123169E-05 == Approx( elements[ t11 ].real() ) );
      CHECK(  8.7747188920057530E-03 == Approx( elements[ t11 ].imag() ) );
      CHECK(  1.7814366984729067E-06 == Approx( elements[ t12 ].real() ) );
      CHECK(  9.0883002678005156E-02 == Approx( elements[ t12 ].imag() ) );
      CHECK( -7.1062891323320817E-04 == Approx( elements[ t13 ].real() ) );
      CHECK(  5.1380747988416891E-05 == Approx( elements[ t13 ].imag() ) );
      CHECK(  1.7814366984729067E-06 == Approx( elements[ t21 ].real() ) );
      CHECK(  9.0883002678005156E-02 == Approx( elements[ t21 ].imag() ) );
      CHECK( -3.0860148985012218E-07 == Approx( elements[ t22 ].real() ) );
      CHECK(  0.94136780169043988 == Approx( elements[ t22 ].imag() ) );
      CHECK(  2.6649964127579322E-05 == Approx( elements[ t23 ].real() ) );
      CHECK(  4.5080252732947884E-04 == Approx( elements[ t23 ].imag() ) );
      CHECK( -7.1062891323320817E-04 == Approx( elements[ t31 ].real() ) );
      CHECK(  5.1380747988416891E-05 == Approx( elements[ t31 ].imag() ) );
      CHECK(  2.6649964127579322E-05 == Approx( elements[ t32 ].real() ) );
      CHECK(  4.5080252732947884E-04 == Approx( elements[ t32 ].imag() ) );
      CHECK( -1.0199940603831528E-02 == Approx( elements[ t33 ].real() ) );
      CHECK(  1.1252686951477420E-04 == Approx( elements[ t33 ].imag() ) );
    } // THEN
  } // GIVEN

  GIVEN( "valid data for a SpinGroup with resonances using negative widths" ) {

    // test based on Pu239 ENDF/B-VIII.0 LRF3 resonance evaluation (all
    // widths used in this test are from the original evaluation)
    // data given in Gamma = 2 gamma^2 P(Er) so conversion is required
    // T matrix values extracted from the samm.f90 routines in NJOY2016.43
    // (note: LRF7 in NJOY2016 doesn't add potential scattering for missing J
    // values)

    // particles
    Particle photon( ParticleID( "g" ), 0.0 * daltons,
                     0.0 * elementary, 1., +1);
    Particle neutron( ParticleID( "n" ), neutronMass,
                      0.0 * elementary, 0.5, +1);
    Particle pu240( ParticleID( "Pu240_e0" ), 2.379916e+2 * neutronMass,
                    94.0 * elementary, 0.5, +1);
    Particle pu239( ParticleID( "Pu239_e0" ), 2.369986e+2 * neutronMass,
                    94.0 * elementary, 0.5, +1);

    // particle pairs
    ParticlePair in( neutron, pu239 );
    ParticlePair out1( photon, pu240 );
    ParticlePair out2( neutron, pu239, ParticlePairID( "fission" ) );

    // channels
    Channel< Photon > capture( in, out1, 0. * electronVolt, { 0, 0.0, 0.0, +1 },
                               { 0.0 * rootBarn },
                               0.0 );
    Channel< Neutron > elastic( in, in, 0. * electronVolt, { 0, 0.5, 0.0, +1 },
                                { 9.410000e-1 * rootBarn },
                                0.0 );
    Channel< Fission > fission1( in, out2, "fission1", 0. * electronVolt,
                                 { 0, 0.0, 0.0, +1 },
                                 { 9.410000e-1 * rootBarn },
                                 0.0 );
    Channel< Fission > fission2( in, out2, "fission2", 0. * electronVolt,
                                 { 0, 0.0, 0.0, +1 },
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
        group( { elastic, fission1, fission2 }, std::move( table ) );

    ReactionChannelID t11( "n,Pu239{0,1/2,0+}->n,Pu239{0,1/2,0+}" );
    ReactionChannelID t12( "n,Pu239{0,1/2,0+}->fission1{0,0,0+}" );
    ReactionChannelID t13( "n,Pu239{0,1/2,0+}->fission2{0,0,0+}" );
    ReactionChannelID t21( "fission1{0,0,0+}->n,Pu239{0,1/2,0+}" );
    ReactionChannelID t22( "fission1{0,0,0+}->fission1{0,0,0+}" );
    ReactionChannelID t23( "fission1{0,0,0+}->fission2{0,0,0+}" );
    ReactionChannelID t31( "fission2{0,0,0+}->n,Pu239{0,1/2,0+}" );
    ReactionChannelID t32( "fission2{0,0,0+}->fission1{0,0,0+}" );
    ReactionChannelID t33( "fission2{0,0,0+}->fission2{0,0,0+}" );

    THEN( "T matrix elements can be calculated" ) {

      std::map< ReactionChannelID, std::complex< double > > elements;
      group.evaluateTMatrix( 1e-5 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK(  8.6094360117330500E-08 == Approx( elements[ t11 ].real() ) );
      CHECK(  1.2969579643023383E-09 == Approx( elements[ t11 ].imag() ) );
      CHECK(  1.2801915235403000E-05 == Approx( elements[ t12 ].real() ) );
      CHECK(  6.5926507156132303E-08 == Approx( elements[ t12 ].imag() ) );
      CHECK(  3.2390394989094051E-05 == Approx( elements[ t13 ].real() ) );
      CHECK(  8.9681494780835694E-07 == Approx( elements[ t13 ].imag() ) );
      CHECK(  1.2801915235403000E-05 == Approx( elements[ t21 ].real() ) );
      CHECK(  6.5926507156132303E-08 == Approx( elements[ t21 ].imag() ) );
      CHECK(  5.9174467948909766E-03 == Approx( elements[ t22 ].real() ) );
      CHECK(  3.7061252598596956E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK( -4.3208179897630039E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( -1.4291551760903298E-05 == Approx( elements[ t23 ].imag() ) );
      CHECK(  3.2390394989094051E-05 == Approx( elements[ t31 ].real() ) );
      CHECK(  8.9681494780835694E-07 == Approx( elements[ t31 ].imag() ) );
      CHECK( -4.3208179897630039E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( -1.4291551760903298E-05 == Approx( elements[ t32 ].imag() ) );
      CHECK(  2.6443852496640871E-02 == Approx( elements[ t33 ].real() ) );
      CHECK(  7.3345304025029671E-04 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e-4 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK(  2.7225547498043351E-07 == Approx( elements[ t11 ].real() ) );
      CHECK(  4.1014378002697183E-09 == Approx( elements[ t11 ].imag() ) );
      CHECK(  2.2765426222498610E-05 == Approx( elements[ t12 ].real() ) );
      CHECK(  1.1724030555007683E-07 == Approx( elements[ t12 ].imag() ) );
      CHECK(  5.7599528678230326E-05 == Approx( elements[ t13 ].real() ) );
      CHECK(  1.5948171792593085E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK(  2.2765426222498610E-05 == Approx( elements[ t21 ].real() ) );
      CHECK(  1.1724030555007683E-07 == Approx( elements[ t21 ].imag() ) );
      CHECK(  5.9174580665652806E-03 == Approx( elements[ t22 ].real() ) );
      CHECK(  3.7061748656442581E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK( -4.3208311458152492E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( -1.4290767857297056E-05 == Approx( elements[ t23 ].imag() ) );
      CHECK(  5.7599528678230326E-05 == Approx( elements[ t31 ].real() ) );
      CHECK(  1.5948171792593085E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK( -4.3208311458152492E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( -1.4290767857297056E-05 == Approx( elements[ t32 ].imag() ) );
      CHECK(  2.6444000478933145E-02 == Approx( elements[ t33 ].real() ) );
      CHECK(  7.3346353581992495E-04 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e-3 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK(  8.6098549373132083E-07 == Approx( elements[ t11 ].real() ) );
      CHECK(  1.2971844538545417E-08 == Approx( elements[ t11 ].imag() ) );
      CHECK(  4.0484071329922845E-05 == Approx( elements[ t12 ].real() ) );
      CHECK(  2.0851550150244872E-07 == Approx( elements[ t12 ].imag() ) );
      CHECK(  1.0243439585180307E-04 == Approx( elements[ t13 ].real() ) );
      CHECK(  2.8364270447639995E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK(  4.0484071329922845E-05 == Approx( elements[ t21 ].real() ) );
      CHECK(  2.0851550150244872E-07 == Approx( elements[ t21 ].imag() ) );
      CHECK(  5.9175708106342824E-03 == Approx( elements[ t22 ].real() ) );
      CHECK(  3.7064286283464446E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK( -4.3209606977718330E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( -1.4289058530520225E-05 == Approx( elements[ t23 ].imag() ) );
      CHECK(  1.0243439585180307E-04 == Approx( elements[ t31 ].real() ) );
      CHECK(  2.8364270447639995E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK( -4.3209606977718330E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( -1.4289058530520225E-05 == Approx( elements[ t32 ].imag() ) );
      CHECK(  2.6445481253788445E-02 == Approx( elements[ t33 ].real() ) );
      CHECK(  7.3355300054415192E-04 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e-2 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK(  2.7238806360778982E-06 == Approx( elements[ t11 ].real() ) );
      CHECK(  4.1071628822804224E-08 == Approx( elements[ t11 ].imag() ) );
      CHECK(  7.2005916425102253E-05 == Approx( elements[ t12 ].real() ) );
      CHECK(  3.7103339178456541E-07 == Approx( elements[ t12 ].imag() ) );
      CHECK(  1.8226982516167094E-04 == Approx( elements[ t13 ].real() ) );
      CHECK(  5.0502838646838534E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK(  7.2005916425102253E-05 == Approx( elements[ t21 ].real() ) );
      CHECK(  3.7103339178456541E-07 == Approx( elements[ t21 ].imag() ) );
      CHECK(  5.9186985672988390E-03 == Approx( elements[ t22 ].real() ) );
      CHECK(  3.7082006184398922E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK( -4.3222503027064754E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( -1.4291349759794570E-05 == Approx( elements[ t23 ].imag() ) );
      CHECK(  1.8226982516167094E-04 == Approx( elements[ t31 ].real() ) );
      CHECK(  5.0502838646838534E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK( -4.3222503027064754E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( -1.4291349759794570E-05 == Approx( elements[ t32 ].imag() ) );
      CHECK(  2.6460301027705024E-02 == Approx( elements[ t33 ].real() ) );
      CHECK(  7.3439943634519798E-04 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e-1 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK(  8.6520079376310224E-06 == Approx( elements[ t11 ].real() ) );
      CHECK(  1.3139947315384008E-07 == Approx( elements[ t11 ].imag() ) );
      CHECK(  1.2829484118662680E-04 == Approx( elements[ t12 ].real() ) );
      CHECK(  6.6234773958743625E-07 == Approx( elements[ t12 ].imag() ) );
      CHECK(  3.2614691456065846E-04 == Approx( elements[ t13 ].real() ) );
      CHECK(  9.0902175262644802E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK(  1.2829484118662680E-04 == Approx( elements[ t21 ].real() ) );
      CHECK(  6.6234773958743625E-07 == Approx( elements[ t21 ].imag() ) );
      CHECK(  5.9300001409376933E-03 == Approx( elements[ t22 ].real() ) );
      CHECK(  3.7235471877830503E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK( -4.3351712695469916E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( -1.4376112601285966E-05 == Approx( elements[ t23 ].imag() ) );
      CHECK(  3.2614691456065846E-04 == Approx( elements[ t31 ].real() ) );
      CHECK(  9.0902175262644802E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK( -4.3351712695469916E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( -1.4376112601285966E-05 == Approx( elements[ t32 ].imag() ) );
      CHECK(  2.6609444431195336E-02 == Approx( elements[ t33 ].real() ) );
      CHECK(  7.4278767834262715E-04 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e+0 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK(  2.8646113295746252E-05 == Approx( elements[ t11 ].real() ) );
      CHECK(  4.6744028655544269E-07 == Approx( elements[ t11 ].imag() ) );
      CHECK(  2.3265382264289994E-04 == Approx( elements[ t12 ].real() ) );
      CHECK(  1.2141299452405830E-06 == Approx( elements[ t12 ].imag() ) );
      CHECK(  6.1843879931493977E-04 == Approx( elements[ t13 ].real() ) );
      CHECK(  1.8291330337424256E-05 == Approx( elements[ t13 ].imag() ) );
      CHECK(  2.3265382264289994E-04 == Approx( elements[ t21 ].real() ) );
      CHECK(  1.2141299452405830E-06 == Approx( elements[ t21 ].imag() ) );
      CHECK(  6.0454441690149666E-03 == Approx( elements[ t22 ].real() ) );
      CHECK(  3.8741401117764872E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK( -4.4689947772132772E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( -1.5494916692727437E-05 == Approx( elements[ t23 ].imag() ) );
      CHECK(  6.1843879931493977E-04 == Approx( elements[ t31 ].real() ) );
      CHECK(  1.8291330337424256E-05 == Approx( elements[ t31 ].imag() ) );
      CHECK( -4.4689947772132772E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( -1.5494916692727437E-05 == Approx( elements[ t32 ].imag() ) );
      CHECK(  2.8201028109277028E-02 == Approx( elements[ t33 ].real() ) );
      CHECK(  8.3472384861434286E-04 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e+1 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK(  1.9449625356215055E-04 == Approx( elements[ t11 ].real() ) );
      CHECK(  1.0372444594343632E-05 == Approx( elements[ t11 ].imag() ) );
      CHECK(  5.1510765829269141E-04 == Approx( elements[ t12 ].real() ) );
      CHECK(  2.1290330431659829E-06 == Approx( elements[ t12 ].imag() ) );
      CHECK(  3.0696316664288220E-03 == Approx( elements[ t13 ].real() ) );
      CHECK(  2.3494912154491140E-04 == Approx( elements[ t13 ].imag() ) );
      CHECK(  5.1510765829269141E-04 == Approx( elements[ t21 ].real() ) );
      CHECK(  2.1290330431659829E-06 == Approx( elements[ t21 ].imag() ) );
      CHECK(  7.5085394656881450E-03 == Approx( elements[ t22 ].real() ) );
      CHECK(  6.0093079962443753E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK( -6.6289359270729114E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( -5.2355066387279181E-05 == Approx( elements[ t23 ].imag() ) );
      CHECK(  3.0696316664288220E-03 == Approx( elements[ t31 ].real() ) );
      CHECK(  2.3494912154491140E-04 == Approx( elements[ t31 ].imag() ) );
      CHECK( -6.6289359270729114E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( -5.2355066387279181E-05 == Approx( elements[ t32 ].imag() ) );
      CHECK(  7.2132075903575107E-02 == Approx( elements[ t33 ].real() ) );
      CHECK(  5.5052043773717698E-03 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e+2 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -1.1372032929517255E-04 == Approx( elements[ t11 ].real() ) );
      CHECK(  5.1934158312836712E-07 == Approx( elements[ t11 ].imag() ) );
      CHECK( -6.3526162811888636E-04 == Approx( elements[ t12 ].real() ) );
      CHECK(  3.5812399114220898E-06 == Approx( elements[ t12 ].imag() ) );
      CHECK( -2.6816325859419470E-04 == Approx( elements[ t13 ].real() ) );
      CHECK(  1.4113606496224005E-06 == Approx( elements[ t13 ].imag() ) );
      CHECK( -6.3526162811888636E-04 == Approx( elements[ t21 ].real() ) );
      CHECK(  3.5812399114220898E-06 == Approx( elements[ t21 ].imag() ) );
      CHECK( -5.3265203913305265E-03 == Approx( elements[ t22 ].real() ) );
      CHECK(  3.0308858408396067E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK(  1.9907302168212490E-04 == Approx( elements[ t23 ].real() ) );
      CHECK( -2.0277682885686721E-06 == Approx( elements[ t23 ].imag() ) );
      CHECK( -2.6816325859419470E-04 == Approx( elements[ t31 ].real() ) );
      CHECK(  1.4113606496224005E-06 == Approx( elements[ t31 ].imag() ) );
      CHECK(  1.9907302168212490E-04 == Approx( elements[ t32 ].real() ) );
      CHECK( -2.0277682885686721E-06 == Approx( elements[ t32 ].imag() ) );
      CHECK( -5.4078837827793739E-03 == Approx( elements[ t33 ].real() ) );
      CHECK(  3.0719306457874704E-05 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1e+3 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -2.3356639937483146E-05 == Approx( elements[ t11 ].real() ) );
      CHECK(  6.8316924521103610E-09 == Approx( elements[ t11 ].imag() ) );
      CHECK( -6.2756282536005842E-05 == Approx( elements[ t12 ].real() ) );
      CHECK(  2.0285942144168103E-08 == Approx( elements[ t12 ].imag() ) );
      CHECK( -4.3916352638395594E-05 == Approx( elements[ t13 ].real() ) );
      CHECK(  2.0780444778526720E-08 == Approx( elements[ t13 ].imag() ) );
      CHECK( -6.2756282536005842E-05 == Approx( elements[ t21 ].real() ) );
      CHECK(  2.0285942144168103E-08 == Approx( elements[ t21 ].imag() ) );
      CHECK( -2.9399211546199931E-04 == Approx( elements[ t22 ].real() ) );
      CHECK(  9.5117202175650674E-08 == Approx( elements[ t22 ].imag() ) );
      CHECK(  1.3971269282865114E-05 == Approx( elements[ t23 ].real() ) );
      CHECK( -7.9320697666663586E-09 == Approx( elements[ t23 ].imag() ) );
      CHECK( -4.3916352638395594E-05 == Approx( elements[ t31 ].real() ) );
      CHECK(  2.0780444778526720E-08 == Approx( elements[ t31 ].imag() ) );
      CHECK(  1.3971269282865114E-05 == Approx( elements[ t32 ].real() ) );
      CHECK( -7.9320697666663586E-09 == Approx( elements[ t32 ].imag() ) );
      CHECK( -4.4949738781199479E-04 == Approx( elements[ t33 ].real() ) );
      CHECK(  2.1349394311927584E-07 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 2e+3 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -1.6226096271959019E-05 == Approx( elements[ t11 ].real() ) );
      CHECK(  2.4077043297474406E-09 == Approx( elements[ t11 ].imag() ) );
      CHECK( -3.6418204959627669E-05 == Approx( elements[ t12 ].real() ) );
      CHECK(  5.9144855333236386E-09 == Approx( elements[ t12 ].imag() ) );
      CHECK( -2.5976368082316346E-05 == Approx( elements[ t13 ].real() ) );
      CHECK(  6.2188579236123019E-09 == Approx( elements[ t13 ].imag() ) );
      CHECK( -3.6418204959627669E-05 == Approx( elements[ t21 ].real() ) );
      CHECK(  5.9144855333236386E-09 == Approx( elements[ t21 ].imag() ) );
      CHECK( -1.4342807020535593E-04 == Approx( elements[ t22 ].real() ) );
      CHECK(  2.3028683929806548E-08 == Approx( elements[ t22 ].imag() ) );
      CHECK(  6.8723149961035684E-06 == Approx( elements[ t23 ].real() ) );
      CHECK( -1.6429830151187082E-09 == Approx( elements[ t23 ].imag() ) );
      CHECK( -2.5976368082316346E-05 == Approx( elements[ t31 ].real() ) );
      CHECK(  6.2188579236123019E-09 == Approx( elements[ t31 ].imag() ) );
      CHECK(  6.8723149961035684E-06 == Approx( elements[ t32 ].real() ) );
      CHECK( -1.6429830151187082E-09 == Approx( elements[ t32 ].imag() ) );
      CHECK( -2.2271662131030396E-04 == Approx( elements[ t33 ].real() ) );
      CHECK(  5.2613057057688388E-08 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 1.541700e+1 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK(  9.8098208745806199E-05 == Approx( elements[ t11 ].real() ) );
      CHECK(  2.5786090151820738E-03 == Approx( elements[ t11 ].imag() ) );
      CHECK(  7.1794528627063480E-04 == Approx( elements[ t12 ].real() ) );
      CHECK( -5.2664390914965053E-05 == Approx( elements[ t12 ].imag() ) );
      CHECK( -2.6735454741440894E-05 == Approx( elements[ t13 ].real() ) );
      CHECK(  4.9399281630055664E-02 == Approx( elements[ t13 ].imag() ) );
      CHECK(  7.1794528627063480E-04 == Approx( elements[ t21 ].real() ) );
      CHECK( -5.2664390914965053E-05 == Approx( elements[ t21 ].imag() ) );
      CHECK(  8.7907704686553713E-03 == Approx( elements[ t22 ].real() ) );
      CHECK(  8.3343078573257687E-05 == Approx( elements[ t22 ].imag() ) );
      CHECK( -6.4495186169242103E-05 == Approx( elements[ t23 ].real() ) );
      CHECK( -1.1400747932943543E-03 == Approx( elements[ t23 ].imag() ) );
      CHECK( -2.6735454741440894E-05 == Approx( elements[ t31 ].real() ) );
      CHECK(  4.9399281630055664E-02 == Approx( elements[ t31 ].imag() ) );
      CHECK( -6.4495186169242103E-05 == Approx( elements[ t32 ].real() ) );
      CHECK( -1.1400747932943543E-03 == Approx( elements[ t32 ].imag() ) );
      CHECK(  1.2151850651597894E-05 == Approx( elements[ t33 ].real() ) );
      CHECK(  0.94658992569665734 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 3.232700e+1 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -2.3172581926204482E-04 == Approx( elements[ t11 ].real() ) );
      CHECK(  4.9392463116775245E-03 == Approx( elements[ t11 ].imag() ) );
      CHECK(  9.3031529979401310E-04 == Approx( elements[ t12 ].real() ) );
      CHECK(  1.2142454746259005E-02 == Approx( elements[ t12 ].imag() ) );
      CHECK( -5.4616126769468443E-04 == Approx( elements[ t13 ].real() ) );
      CHECK( -5.9913588141388432E-02 == Approx( elements[ t13 ].imag() ) );
      CHECK(  9.3031529979401310E-04 == Approx( elements[ t21 ].real() ) );
      CHECK(  1.2142454746259005E-02 == Approx( elements[ t21 ].imag() ) );
      CHECK(  1.6657407521629055E-02 == Approx( elements[ t22 ].real() ) );
      CHECK(  3.0080787789397462E-02 == Approx( elements[ t22 ].imag() ) );
      CHECK(  1.7964258198506294E-03 == Approx( elements[ t23 ].real() ) );
      CHECK( -0.14712700890740524 == Approx( elements[ t23 ].imag() ) );
      CHECK( -5.4616126769468443E-04 == Approx( elements[ t31 ].real() ) );
      CHECK( -5.9913588141388432E-02 == Approx( elements[ t31 ].imag() ) );
      CHECK(  1.7964258198506294E-03 == Approx( elements[ t32 ].real() ) );
      CHECK( -0.14712700890740524 == Approx( elements[ t32 ].imag() ) );
      CHECK( -1.2795632018180753E-03 == Approx( elements[ t33 ].real() ) );
      CHECK(  0.72742940272800227 == Approx( elements[ t33 ].imag() ) );

      group.evaluateTMatrix( 4.753400e+1 * electronVolt, elements );
      CHECK( 9 == elements.size() );
      CHECK( -7.6957724678165175E-05 == Approx( elements[ t11 ].real() ) );
      CHECK(  8.7745001866921350E-03 == Approx( elements[ t11 ].imag() ) );
      CHECK(  3.8148891044013242E-06 == Approx( elements[ t12 ].real() ) );
      CHECK(  9.0882977365792011E-02 == Approx( elements[ t12 ].imag() ) );
      CHECK( -5.0676531201553228E-04 == Approx( elements[ t13 ].real() ) );
      CHECK( -3.5327861901406070E-05 == Approx( elements[ t13 ].imag() ) );
      CHECK(  3.8148891044013242E-06 == Approx( elements[ t21 ].real() ) );
      CHECK(  9.0882977365792011E-02 == Approx( elements[ t21 ].imag() ) );
      CHECK( -4.5182942087391158E-07 == Approx( elements[ t22 ].real() ) );
      CHECK(  0.94136780800441566 == Approx( elements[ t22 ].imag() ) );
      CHECK(  8.2692928352256600E-05 == Approx( elements[ t23 ].real() ) );
      CHECK( -4.5245392366873201E-04 == Approx( elements[ t23 ].imag() ) );
      CHECK( -5.0676531201553228E-04 == Approx( elements[ t31 ].real() ) );
      CHECK( -3.5327861901406070E-05 == Approx( elements[ t31 ].imag() ) );
      CHECK(  8.2692928352256600E-05 == Approx( elements[ t32 ].real() ) );
      CHECK( -4.5245392366873201E-04 == Approx( elements[ t32 ].imag() ) );
      CHECK( -1.5953943874582942E-02 == Approx( elements[ t33 ].real() ) );
      CHECK(  2.6830949543030767E-04 == Approx( elements[ t33 ].imag() ) );
    } // THEN
  } // GIVEN

  GIVEN( "valid data for a SpinGroup with one eliminated capture channel, "
         "one elastic channel and a proton channel" ) {

    // test based on Cl35 ENDF/B-VIII.0 LRF7 resonance evaluation
    // data given in Gamma = 2 gamma^2 P(Er) so conversion is required
    // cross section values extracted from NJOY2016.43

    // because the oribital angular momentum l = 0 for these SpinGroup,
    // SpinGroup< ReichMoore, ShiftFactor > and SpinGroup< ReichMoore, Constant >
    // should give the same results

    // using SpinGroup< ReichMoore, ShiftFactor > is equivalent to NJOY2016's
    // LRF7 reconstruction

    // particles
    Particle photon( ParticleID( "g" ), 0.0 * daltons,
                     0.0 * elementary, 1., +1);
    Particle neutron( ParticleID( "n" ), neutronMass,
                      0.0 * elementary, 0.5, +1);
    Particle proton( ParticleID( "p" ), 9.986235e-1 * neutronMass,
                     elementary, 0.5, +1);
    Particle cl36( ParticleID( "Cl36_e0" ), 3.565932e+1 * neutronMass,
                   17.0 * elementary, 0., +1);
    Particle cl35( ParticleID( "Cl35_e0" ), 3.466845e+1 * neutronMass,
                   17.0 * elementary, 1.5, +1);
    Particle s36( ParticleID( "S36_e0" ), 3.466863e+1 * neutronMass,
                  16.0 * elementary, 1.5, +1);

    // particle pairs
    ParticlePair in( neutron, cl35 );
    ParticlePair out1( photon, cl36 );
    ParticlePair out2( proton, s36 );

    // channels
    Channel< Photon > capture( in, out1, 0. * electronVolt, { 0, 0.0, 1.0, +1 },
                               { 0.0 * rootBarn },
                               0.0 );
    Channel< Neutron > elastic( in, in, 0. * electronVolt, { 0, 1.0, 1.0, +1 },
                                { 4.822220e-1 * rootBarn,
                                  3.667980e-1 * rootBarn },
                                0.0 );
    Channel< ChargedParticle > protonemission( in, out2,
                                               6.152200e+5 * electronVolt,
                                               { 0, 1.0, 1.0, +1 },
                                               { 4.822220e-1 * rootBarn,
                                                 3.667980e-1 * rootBarn },
                                               0.0 );

    // conversion from Gamma to gamma
    auto eGamma = [&] ( double width, const Energy& energy ) -> ReducedWidth {
      return std::sqrt( width / 2. / elastic.penetrability( energy ) ) *
             rootElectronVolt;
    };
    auto pGamma = [&] ( double width, const Energy& energy ) -> ReducedWidth {
      return std::sqrt( width / 2. / protonemission.penetrability( energy ) ) *
             rootElectronVolt;
    };
    auto cGamma = [&] ( double width ) -> ReducedWidth {
      return std::sqrt( width / 2. ) * rootElectronVolt;
    };

    // single resonance table
    ResonanceTable single(
      { elastic.channelID(), protonemission.channelID() },
      { Resonance( 6.823616e+4 * electronVolt,
                   { eGamma( 2.179040e+2, 6.823616e+4 * electronVolt ),
                     pGamma( 1.000000e-5, 6.823616e+4 * electronVolt ) },
                   cGamma( 3.933600e-1 ) ) } );
    ResonanceTable single2 = single;

    // multiple resonance table
    ResonanceTable multiple(
      { elastic.channelID(), protonemission.channelID() },
      { Resonance( 6.823616e+4 * electronVolt,
                   { eGamma( 2.179040e+2, 6.823616e+4 * electronVolt ),
                     pGamma( 1.000000e-5, 6.823616e+4 * electronVolt ) },
                   cGamma( 3.933600e-1 ) ),
        Resonance( 1.825230e+5 * electronVolt,
                   { eGamma( 1.759740e+3, 1.825230e+5 * electronVolt ),
                     pGamma( 4.000000e-1, 1.825230e+5 * electronVolt ) },
                   cGamma( 7.451500e-1 ) ),
        Resonance( 2.397427e+5 * electronVolt,
                   { eGamma( 2.685470e+2, 2.397427e+5 * electronVolt ),
                     pGamma( 0.0, 2.397427e+5 * electronVolt ) },
                   cGamma( 6.871600e-1 ) ) } );
    ResonanceTable multiple2 = multiple;

    SpinGroup< ReichMoore, ShiftFactor >
        group1( { elastic, protonemission }, std::move( single ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group2( { elastic, protonemission }, std::move( multiple ) );
    SpinGroup< ReichMoore, Constant >
        group3( { elastic, protonemission }, std::move( single2 ) );
    SpinGroup< ReichMoore, Constant >
        group4( { elastic, protonemission }, std::move( multiple2 ) );

    ReactionChannelID t11( "n,Cl35{0,1,1+}->n,Cl35{0,1,1+}" );
    ReactionChannelID t12( "n,Cl35{0,1,1+}->p,S36{0,1,1+}" );
    ReactionChannelID t21( "p,S36{0,1,1+}->n,Cl35{0,1,1+}" );
    ReactionChannelID t22( "p,S36{0,1,1+}->p,S36{0,1,1+}" );

    THEN( "cross sections can be calculated for a single resonance using the "
          "ShiftFactor boundary condition" ) {

      std::map< ReactionChannelID, std::complex< double > > elements;
      group1.evaluateTMatrix( 1e-5 * electronVolt, elements );
      CHECK( 4 == elements.size() );
      CHECK(  1.9329175916735239E-08 == Approx( elements[ t11 ].real() ) );
      CHECK(  5.6087466302805892E-14 == Approx( elements[ t11 ].imag() ) );
      CHECK(  7.3449798785626636E-10 == Approx( elements[ t12 ].real() ) );
      CHECK(  2.1312926800826200E-15 == Approx( elements[ t12 ].imag() ) );
      CHECK(  7.3449798785626636E-10 == Approx( elements[ t21 ].real() ) );
      CHECK(  2.1312926800826200E-15 == Approx( elements[ t21 ].imag() ) );
      CHECK(  2.7910517059230386E-11 == Approx( elements[ t22 ].real() ) );
      CHECK(  8.0987942362204981E-17 == Approx( elements[ t22 ].imag() ) );
    } // THEN

    THEN( "cross sections can be calculated for multiple resonances using the "
          "ShiftFactor boundary condition" ) {

      std::map< ReactionChannelID, std::complex< double > > elements;
      group2.evaluateTMatrix( 1e-5 * electronVolt, elements );
      CHECK( 4 == elements.size() );
      CHECK(  5.8627792036991544E-08 == Approx( elements[ t11 ].real() ) );
      CHECK(  1.4120321860054129E-13 == Approx( elements[ t11 ].imag() ) );
      CHECK(  6.3514819342570394E-08 == Approx( elements[ t12 ].real() ) );
      CHECK(  1.4100873250164967E-13 == Approx( elements[ t12 ].imag() ) );
      CHECK(  6.3514819342570394E-08 == Approx( elements[ t21 ].real() ) );
      CHECK(  1.4100873250164967E-13 == Approx( elements[ t21 ].imag() ) );
      CHECK(  1.1048788944230737E-07 == Approx( elements[ t22 ].real() ) );
      CHECK(  2.4179854102688846E-13 == Approx( elements[ t22 ].imag() ) );
    } // THEN

    THEN( "cross sections can be calculated for a single resonance using the "
          "Constant boundary condition" ) {

      std::map< ReactionChannelID, std::complex< double > > elements;
      group3.evaluateTMatrix( 1e-5 * electronVolt, elements );
      CHECK( 4 == elements.size() );
      CHECK(  1.9329175916735239E-08 == Approx( elements[ t11 ].real() ) );
      CHECK(  5.6087466302805892E-14 == Approx( elements[ t11 ].imag() ) );
      CHECK(  7.3449798785626636E-10 == Approx( elements[ t12 ].real() ) );
      CHECK(  2.1312926800826200E-15 == Approx( elements[ t12 ].imag() ) );
      CHECK(  7.3449798785626636E-10 == Approx( elements[ t21 ].real() ) );
      CHECK(  2.1312926800826200E-15 == Approx( elements[ t21 ].imag() ) );
      CHECK(  2.7910517059230386E-11 == Approx( elements[ t22 ].real() ) );
      CHECK(  8.0987942362204981E-17 == Approx( elements[ t22 ].imag() ) );
    } // THEN

    THEN( "cross sections can be calculated for multiple resonances using the "
          "Constant boundary condition" ) {

      std::map< ReactionChannelID, std::complex< double > > elements;
      group4.evaluateTMatrix( 1e-5 * electronVolt, elements );
      CHECK( 4 == elements.size() );
      CHECK(  5.8627792036991544E-08 == Approx( elements[ t11 ].real() ) );
      CHECK(  1.4120321860054129E-13 == Approx( elements[ t11 ].imag() ) );
      CHECK(  6.3514819342570394E-08 == Approx( elements[ t12 ].real() ) );
      CHECK(  1.4100873250164967E-13 == Approx( elements[ t12 ].imag() ) );
      CHECK(  6.3514819342570394E-08 == Approx( elements[ t21 ].real() ) );
      CHECK(  1.4100873250164967E-13 == Approx( elements[ t21 ].imag() ) );
      CHECK(  1.1048788944230737E-07 == Approx( elements[ t22 ].real() ) );
      CHECK(  2.4179854102688846E-13 == Approx( elements[ t22 ].imag() ) );
    } // THEN
  } // GIVEN
} // SCENARIO
