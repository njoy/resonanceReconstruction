#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;
using namespace dimwits;

breitWigner::lvalue::Type makeLvalue( int l );

namespace {

struct SLBW : breitWigner::singleLevel::Base {
  using breitWigner::singleLevel::Base::psiChi;
};

};

SCENARIO("call"){
  auto base = makeLvalue( 0 );

  const auto& lValue =
    static_cast< const breitWigner::singleLevel::Lvalue& >( base );
  
  auto awr = 2.360045E+2;
  
  auto a = channelRadius( awr );
  auto ap = radius( 9.309000E-1 );
  auto k = neutronWaveNumber( awr );
  auto rho = [&]( auto energy ) -> double { return a(energy) * k(energy); };

  auto energies =
    ranges::view::linear_distribute( 1.0E-5, 200., 100000 )
    | ranges::view::transform( [&]( double d ){ return d * electronVolts; } );

  auto evaluate = [&]( auto&& energy ){
    auto psiChi = SLBW::psiChi( energy );
    auto channelRatio = a(energy) * k(energy);
    auto scatteringRatio = ap(energy) * k(energy);
    return std::get<2>( lValue( energy,
                                psiChi,
                                channelRatio,
                                scatteringRatio,
                                rho ).data );
  };

  auto unwrap = []( auto q ){ return q.value; };
  std::cout << ( energies | ranges::view::transform( unwrap ) ) << std::endl;

  auto crossSection = energies | ranges::view::transform( evaluate );
  std::cout << crossSection << std::endl;

  /* this isn't really a test of correctness,
   * but does establish that the template instantiations succeed
   */
}
