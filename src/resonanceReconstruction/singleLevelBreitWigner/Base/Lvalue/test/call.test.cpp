#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction::singleLevelBreitWigner;
using namespace dimwits;

Base::Lvalue makeLvalue();

SCENARIO("call"){
  auto lValue = makeLvalue();
  auto awr = 2.360045E+2;
  auto l = 0;
  
  auto a = channelRadius( awr );
  auto ap = radius( 9.309000E-1 );
  auto k = neutronWaveNumber( awr );
  auto phi = [&]( auto energy ){ return a(energy) * k(energy); };

  auto energies =
    ranges::view::linear_distribute( 1.0E-5, 200., 100000 )
    | ranges::view::transform( [&]( double d ){ return d * electronVolts; } );

  auto evaluate = [&]( auto&& energy ){
    auto psiChi = Base::psiChi( energy );
    auto channelRatio = a(energy) * k(energy);
    auto scatteringRatio = ap(energy) * k(energy);
    return std::get<2>( lValue( energy,
                                psiChi,
                                channelRatio,
                                scatteringRatio,
                                a ).data );
  };

  auto unwrap = []( auto q ){ return q.value; };
  std::cout << ( energies | ranges::view::transform( unwrap ) ) << std::endl;

  auto crossSection = energies | ranges::view::transform( evaluate );
  std::cout << crossSection << std::endl;

  /* this isn't really a test of correctness,
   * but does establish that the template instantiations succeed
   */
}
