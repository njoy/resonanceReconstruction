#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction/rmatrix/LMatrixCalculator.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
template < typename BoundaryOption > using LMatrixCalculator =
    rmatrix::LMatrixCalculator< BoundaryOption >;
using Constant = rmatrix::Constant;

using Particle = rmatrix::Particle;
using ParticlePair = rmatrix::ParticlePair;
using ParticleID = rmatrix::ParticleID;
using ParticlePairID = rmatrix::ParticlePairID;
template < typename Type > using Channel = rmatrix::Channel< Type >;
using ChannelID = rmatrix::ChannelID;
using ChannelQuantumNumbers = rmatrix::ChannelQuantumNumbers;
using ChannelRadii = rmatrix::ChannelRadii;
using BoundaryCondition = rmatrix::BoundaryCondition;

constexpr ElectricalCharge elementary = dimwits::constant::elementaryCharge;

struct Test1 {};
struct Test2 {};
struct Test3 {};

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {

template <> inline
double calculatePenetrability< Test1 >(
  const unsigned int, const double, const double ) { return 1.; }

template <> inline
double calculatePenetrability< Test2 >(
  const unsigned int, const double, const double ) { return 2.; }

template <> inline
double calculatePenetrability< Test3 >(
  const unsigned int, const double, const double ) { return 3.; }

template <> inline
double calculateShiftFactor< Test1 >(
  const unsigned int, const double, const double ) { return 3.; }

template <> inline
double calculateShiftFactor< Test2 >(
  const unsigned int, const double, const double ) { return 2.; }

template <> inline
double calculateShiftFactor< Test3 >(
  const unsigned int, const double, const double ) { return 1.; }
}
}
}

using TestChannel = std::variant< Channel< Test1 >,
                                  Channel< Test2 >,
                                  Channel< Test3 > >;

std::vector< TestChannel > createChannels();

SCENARIO( "LMatrixCalculator< Constant >" ) {

  GIVEN( "valid data for channels" ) {

    Energy energy = 1e-5 * electronVolts;

    const auto channels = createChannels();
    const auto penetrabilities = channels
        | ranges::views::transform(
            [&] ( const auto& channel )
                { return std::visit(
                         [&] ( const auto& channel )
                             { return channel.penetrability( energy ); },
                         channel ); } );

    THEN( "a calculator can be constructed" ) {

      LMatrixCalculator< Constant > calculator( channels.size() );

      CHECK( 3 == calculator.matrix().diagonal().size() );
      CHECK( 0.0 == calculator.matrix().diagonal()[0] );
      CHECK( 0.0 == calculator.matrix().diagonal()[1] );
      CHECK( 0.0 == calculator.matrix().diagonal()[2] );
    } // THEN

    THEN( "L = S - B + iP can be calculated" ) {

      LMatrixCalculator< Constant > calculator( channels.size() );

      decltype(auto) matrix = calculator( energy, penetrabilities, channels );

      CHECK( 3 == matrix.diagonal().size() );
      CHECK( std::complex< double >( 3., 1. ) == matrix.diagonal()[0] );
      CHECK( std::complex< double >( 2., 2. ) == matrix.diagonal()[1] );
      CHECK( std::complex< double >( 1., 3. ) == matrix.diagonal()[2] );
    } // THEN
  } // GIVEN
} // SCENARIO

std::vector< TestChannel > createChannels() {

  // particles
  Particle photon( ParticleID( "g" ), 0.0 * daltons, 0.0 * coulombs, 1., +1);
  Particle neutron( ParticleID( "n" ), 1.00866491582 * daltons,
                    0.0 * coulombs, 0.5, +1);
  Particle proton( ParticleID( "p" ), 1.00727647 * daltons,
                   elementary, 0.5, +1);
  Particle cl36( ParticleID( "Cl36_e0" ), 35.968306822 * daltons,
                 17.0 * elementary, 0., +1);
  Particle cl35( ParticleID( "Cl35_e0" ), 34.968852694 * daltons,
                 17.0 * elementary, 1.5, +1);
  Particle cl35_e1( ParticleID( "Cl35_e1" ), 34.968852694 * daltons,
                    17.0 * elementary, 1.5, +1);
  Particle s36( ParticleID( "S36_e0" ), 35.967080699 * daltons,
                16.0 * elementary, 1.5, +1);

  // particle pairs
  ParticlePair elasticPair( neutron, cl35 );
  ParticlePair inelasticPair( neutron, cl35_e1 );
  ParticlePair capturePair( photon, cl36 );
  ParticlePair protonEmissionPair( proton, s36 );

  // Q values
  QValue elasticQ = 0.0 * electronVolt;
  QValue inelasticQ = -1.219440e+6 * electronVolt;
  QValue captureQ = 0.0 * electronVolt;
  QValue protonEmissionQ = 6.152200e+5 * electronVolt;

  // quantum numbers
  ChannelQuantumNumbers elasticNumbers( 0, 1.0, 1.0, +1 );
  ChannelQuantumNumbers inelasticNumbers( 0, 1.0, 1.0, +1 );
  ChannelQuantumNumbers captureNumbers( 0, 0.0, 1.0, +1 );
  ChannelQuantumNumbers protonEmissionNumbers( 0, 1.0, 1.0, +1 );

  // channel radii
  ChannelRadii elasticRadii( 4.822220e-1 * rootBarn,
                             3.667980e-1 * rootBarn );
  ChannelRadii inelasticRadii( 4.822220e-1 * rootBarn,
                             3.667980e-1 * rootBarn );
  ChannelRadii captureRadii( 0.0 * rootBarn );
  ChannelRadii protonEmissionRadii( 4.822220e-1 * rootBarn,
                                    3.667980e-1 * rootBarn );

  return { Channel< Test1 >( elasticPair, elasticPair, elasticQ,
                             elasticNumbers, elasticRadii ),
           Channel< Test2 >( inelasticPair, inelasticPair, inelasticQ,
                             inelasticNumbers, inelasticRadii ),
           Channel< Test3 >( protonEmissionPair, protonEmissionPair,
                             protonEmissionQ,
                             protonEmissionNumbers, protonEmissionRadii ) };
}
