#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction/rmatrix/RLMatrixCalculator.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
template < typename Formalism, typename BoundaryOption > using RLMatrixCalculator =
    rmatrix::RLMatrixCalculator< Formalism, BoundaryOption >;
using Constant = rmatrix::Constant;

using Particle = rmatrix::Particle;
using ParticlePair = rmatrix::ParticlePair;
using ParticleID = rmatrix::ParticleID;
using ParticlePairID = rmatrix::ParticlePairID;
using Neutron = rmatrix::Neutron;
using Photon = rmatrix::Photon;
using Fission = rmatrix::Fission;
using ChargedParticle = rmatrix::ChargedParticle;
template < typename Type > using Channel = rmatrix::Channel< Type >;
using ParticleChannel = rmatrix::ParticleChannel;
using ChannelID = rmatrix::ChannelID;
using ChannelQuantumNumbers = rmatrix::ChannelQuantumNumbers;
using ChannelRadii = rmatrix::ChannelRadii;
using BoundaryCondition = rmatrix::BoundaryCondition;
using Resonance = rmatrix::Resonance;
using ResonanceTable = rmatrix::ResonanceTable;
using ShiftFactor = rmatrix::ShiftFactor;
using Constant = rmatrix::Constant;
using ReichMoore = rmatrix::ReichMoore;

constexpr AtomicMass neutronMass = 1.008664 * daltons;
constexpr ElectricalCharge elementary = dimwits::constant::elementaryCharge;
constexpr double e = 1.6021766208e-19;

SCENARIO( "evaluate" ) {

  GIVEN( "valid channel and resonance table data" ) {

    // test based on Cl35 ENDF/B-VIII.0 LRF7 resonance evaluation
    // data given in Gamma = 2 gamma^2 P(Er) so conversion is required

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
    Particle cl36( ParticleID( "Cl36" ), 3.565932e+1 * neutronMass,
                   17.0 * elementary, 0., +1);
    Particle cl35( ParticleID( "Cl35" ), 3.466845e+1 * neutronMass,
                   17.0 * elementary, 1.5, +1);
    Particle s36( ParticleID( "S36" ), 3.466863e+1 * neutronMass,
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
                                               -1.87263536e+00 );

    // note: the boundary condition for the proton emission channel is set to
    // be equal to the shift factor S at 1e-5 eV. When changing the coulomb
    // wave functions, this value may have to be changed

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

    Energy energy = 1e-5 * electronVolts;

    std::vector< ParticleChannel > channels = { elastic, protonemission };
    const auto penetrabilities = channels
        | ranges::views::transform(
            [&] ( const auto& channel )
                { return std::visit(
                         [&] ( const auto& channel )
                             { return channel.penetrability( energy ); },
                         channel ); } );

 const auto shifts = channels
     | ranges::views::transform(
         [&] ( const auto& channel )
             { return std::visit(
                      [&] ( const auto& channel )
                          { return channel.shiftFactor( energy ); },
                      channel ); } );


    THEN( "a calculator can be constructed" ) {

      RLMatrixCalculator< ReichMoore, ShiftFactor > calculator1( single );
      RLMatrixCalculator< ReichMoore, ShiftFactor > calculator2( multiple );
      RLMatrixCalculator< ReichMoore, Constant > calculator3( single );
      RLMatrixCalculator< ReichMoore, Constant > calculator4( multiple );

      CHECK( 2 == calculator1.matrix().rows() );
      CHECK( 2 == calculator1.matrix().cols() );
      CHECK( 0.0 == calculator1.matrix()(0,0) );
      CHECK( 0.0 == calculator1.matrix()(0,1) );
      CHECK( 0.0 == calculator1.matrix()(1,0) );
      CHECK( 0.0 == calculator1.matrix()(1,1) );

      CHECK( 2 == calculator2.matrix().rows() );
      CHECK( 2 == calculator2.matrix().cols() );
      CHECK( 0.0 == calculator2.matrix()(0,0) );
      CHECK( 0.0 == calculator2.matrix()(0,1) );
      CHECK( 0.0 == calculator2.matrix()(1,0) );
      CHECK( 0.0 == calculator2.matrix()(1,1) );

      CHECK( 2 == calculator3.matrix().rows() );
      CHECK( 2 == calculator3.matrix().cols() );
      CHECK( 0.0 == calculator3.matrix()(0,0) );
      CHECK( 0.0 == calculator3.matrix()(0,1) );
      CHECK( 0.0 == calculator3.matrix()(1,0) );
      CHECK( 0.0 == calculator3.matrix()(1,1) );

      CHECK( 2 == calculator4.matrix().rows() );
      CHECK( 2 == calculator4.matrix().cols() );
      CHECK( 0.0 == calculator4.matrix()(0,0) );
      CHECK( 0.0 == calculator4.matrix()(0,1) );
      CHECK( 0.0 == calculator4.matrix()(1,0) );
      CHECK( 0.0 == calculator4.matrix()(1,1) );
    } // THEN

    THEN( "R_L = ( I - RL )^-1 R can be calculated" ) {

      RLMatrixCalculator< ReichMoore, ShiftFactor > calculator1( single );
      RLMatrixCalculator< ReichMoore, ShiftFactor > calculator2( multiple );
      RLMatrixCalculator< ReichMoore, Constant > calculator3( single );
      RLMatrixCalculator< ReichMoore, Constant > calculator4( multiple );

      decltype(auto) matrix1 = calculator1( energy, single,
                                            penetrabilities, channels );
      decltype(auto) matrix2 = calculator2( energy, multiple,
                                            penetrabilities, channels );
      decltype(auto) matrix3 = calculator3( energy, single,
                                            penetrabilities, channels );
      decltype(auto) matrix4 = calculator4( energy, multiple,
                                            penetrabilities, channels );

      CHECK( 5.93641544e-03 == Approx( matrix1(0,0).real() ) );
      CHECK( 1.72256956e-08 == Approx( matrix1(0,0).imag() ) );
      CHECK( 7.55420060e-05 == Approx( matrix1(0,1).real() ) );
      CHECK( 2.19200225e-10 == Approx( matrix1(0,1).imag() ) );
      CHECK( 7.55420060e-05 == Approx( matrix1(1,0).real() ) );
      CHECK( 2.19200225e-10 == Approx( matrix1(1,0).imag() ) );
      CHECK( 9.61286272e-07 == Approx( matrix1(1,1).real() ) );
      CHECK( 2.78936420e-12 == Approx( matrix1(1,1).imag() ) );

      CHECK( 1.80058856e-02 == Approx( matrix2(0,0).real() ) );
      CHECK( 4.33704410e-08 == Approx( matrix2(0,0).imag() ) );
      CHECK( 6.54248974e-03 == Approx( matrix2(0,1).real() ) );
      CHECK( 1.45270867e-08 == Approx( matrix2(0,1).imag() ) );
      CHECK( 6.54248974e-03 == Approx( matrix2(1,0).real() ) );
      CHECK( 1.45270867e-08 == Approx( matrix2(1,0).imag() ) );
      CHECK( 3.81728883e-03 == Approx( matrix2(1,1).real() ) );
      CHECK( 8.35530059e-09 == Approx( matrix2(1,1).imag() ) );

      CHECK( 5.93641544e-03 == Approx( matrix3(0,0).real() ) );
      CHECK( 1.72256956e-08 == Approx( matrix3(0,0).imag() ) );
      CHECK( 7.55420060e-05 == Approx( matrix3(0,1).real() ) );
      CHECK( 2.19200225e-10 == Approx( matrix3(0,1).imag() ) );
      CHECK( 7.55420060e-05 == Approx( matrix3(1,0).real() ) );
      CHECK( 2.19200225e-10 == Approx( matrix3(1,0).imag() ) );
      CHECK( 9.61286272e-07 == Approx( matrix3(1,1).real() ) );
      CHECK( 2.78936420e-12 == Approx( matrix3(1,1).imag() ) );

      CHECK( 1.80058856e-02 == Approx( matrix4(0,0).real() ) );
      CHECK( 4.33704410e-08 == Approx( matrix4(0,0).imag() ) );
      CHECK( 6.54248974e-03 == Approx( matrix4(0,1).real() ) );
      CHECK( 1.45270867e-08 == Approx( matrix4(0,1).imag() ) );
      CHECK( 6.54248974e-03 == Approx( matrix4(1,0).real() ) );
      CHECK( 1.45270867e-08 == Approx( matrix4(1,0).imag() ) );
      CHECK( 3.81728883e-03 == Approx( matrix4(1,1).real() ) );
      CHECK( 8.35530059e-09 == Approx( matrix4(1,1).imag() ) );
    } // THEN
  } // WHEN
} // SCENARIO
