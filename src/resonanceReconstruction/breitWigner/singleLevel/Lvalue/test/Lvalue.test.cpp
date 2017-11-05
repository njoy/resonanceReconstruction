#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;
using namespace dimwits;

breitWigner::lvalue::Type makeLvalue( int l ){
  auto awr = 2.360045E+2;

  auto a = channelRadius( awr );
  // unused but nice for reference
  // auto ap = radius( 9.309000E-1 );
  // auto k = neutronWaveNumber( awr );

  std::vector< double > resonanceBlock =
  {  -1.000000E+1, 5.000000E-1, 6.217000E-2, 1.581000E-2, 4.500000E-2, 1.360000E-3,
     -4.000000E-1, 5.000000E-1, 4.670000E-2, 3.400000E-4, 4.500000E-2, 1.360000E-3,
      2.855000E+0, 5.000000E-1, 3.808600E-2, 7.470000E-5, 3.680000E-2, 1.211000E-3,
      9.975000E+0, 5.000000E-1, 3.721300E-2, 2.084000E-4, 3.024000E-2, 6.765000E-3,
      1.856000E+1, 5.000000E-1, 4.249000E-2, 3.490000E-3, 3.739000E-2, 1.610000E-3,
      5.980000E+1, 5.000000E-1, 3.930000E-2, 1.550000E-3, 3.480000E-2, 2.950000E-3,
      7.010000E+1, 5.000000E-1, 4.504700E-2, 2.512000E-3, 3.480000E-2, 7.735000E-3,
      8.300000E+1, 5.000000E-1, 5.835600E-2, 1.913000E-2, 3.480000E-2, 4.426000E-3,
      1.100000E+2, 5.000000E-1, 4.883500E-2, 5.245000E-3, 3.480000E-2, 8.790000E-3,
      1.138000E+2, 5.000000E-1, 5.092000E-2, 1.035000E-2, 3.480000E-2, 5.770000E-3,
      1.186000E+2, 5.000000E-1, 6.907000E-2, 3.191000E-2, 3.480000E-2, 2.360000E-3,
      1.225000E+2, 5.000000E-1, 7.021000E-2, 2.667000E-2, 3.480000E-2, 8.740000E-3,
      1.510000E+2, 5.000000E-1, 6.814000E-2, 1.352000E-2, 3.480000E-2, 1.982000E-2,
      1.712000E+2, 5.000000E-1, 1.002200E-1, 6.542000E-2, 3.480000E-2, 0.000000E+0,
      1.825000E+2, 5.000000E-1, 8.359000E-2, 3.242000E-2, 3.480000E-2, 1.637000E-2,
      1.920000E+2, 5.000000E-1, 1.060400E-1, 4.851000E-2, 3.480000E-2, 2.273000E-2 };

  const auto processResonance = [&]( auto&& chunk ){
    const auto energy       = chunk[0] * electronVolts;
    const auto neutronWidth = chunk[3] * electronVolts;
    const auto captureWidth = chunk[4] * electronVolts;
    const auto fissionWidth = chunk[5] * electronVolts;

    const auto waveNumber = neutronWaveNumber( awr );

    const auto weightedCompetitiveWidth = 0.0 * electronVolts;

    const auto channelRadius = a( std::abs(energy) );
    const auto channelRatio  = channelRadius * waveNumber( energy );

    const auto ps = penetrationShift( l, channelRatio );
    const auto& penetrationFactor = ps[0];
    const auto& shiftFactor = ps[1];
    const auto statisticalFactor = ( 2 * chunk[1] + 1 ) / ( 2. );
    return breitWigner::resonance::Type{ energy,
                                         neutronWidth,
                                         captureWidth,
                                         fissionWidth,
                                         weightedCompetitiveWidth,
                                         1. / penetrationFactor,
                                         shiftFactor,
                                         statisticalFactor };
  };

  auto resonances = resonanceBlock
                    | ranges::view::chunk( 6 )
                    | ranges::view::transform( processResonance )
                    | ranges::to_vector;

  return { std::move(resonances), l, 0.0 };
}
