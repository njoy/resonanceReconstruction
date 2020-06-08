#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;
using namespace dimwits;

/*
namespace {

struct SLBW : breitWigner::singleLevel::Base {
  using breitWigner::singleLevel::Base::psiChi;
};

};
*/

SCENARIO("a MLBW resonance reconstruction"){
  auto a = channelRadius( 2.360045E+2 );
  // unused by nice for reference
  // auto ap = radius( 9.309000E-1 );
  auto k = neutronWaveNumber( 2.360045E+2 );
  int l = 0;
  
  auto evaluate = [&]( auto&& resonance, auto energy ){
    auto waveNumber = k(energy);
    auto channelRatio = a(energy) * waveNumber;
    auto scatteringRatio = a(energy) * waveNumber;
    const auto ps = penetrationShift( l, channelRatio );
    const auto& penetrationFactor = ps[0];
    const auto& shiftFactor = ps[1];
    
    const auto phase = phaseShift( l, scatteringRatio );
    
    const auto trig = [&]() -> std::array<double, 3> {
      const auto sin = std::sin( phase );
      const auto sinSquared = sin * sin;
      const auto sin2 = std::sin( 2. * phase );
      const auto cos2 = std::cos( 2. * phase );
      return {{ sinSquared, sin2, cos2 }};
    }();

    auto psiChi = breitWigner::psiChi( energy );

    return resonance( penetrationFactor,
                      shiftFactor,
                      trig[0],
                      trig[1],
                      0.0 * electronVolts,
                      psiChi );
  };
  
  SECTION( "third resonance of pu238 of ENDF 6" ){
    auto resonanceEnergy = 2.885 * electronVolts;
    auto totalWidth = 3.8086E-2 * electronVolts;
    auto neutronWidth = 7.47E-5 * electronVolts;
    auto captureWidth = 3.68E-2 * electronVolts;
    auto fissionWidth = 1.211E-3 * electronVolts;
    auto weightedCompetitiveWidth = 0. * electronVolts;

    auto ps = penetrationShift( l,
                                k( resonanceEnergy )
                                * a( std::abs( resonanceEnergy ) ) );
  
    auto inversePenetrationFactor = 1. / ps[0];
    auto shiftFactor = ps[1];
  
    auto j = 0.5;
    auto i = 0.0;
    auto statisticalFactor = ( 2. * j + 1. ) / ( 4. * i + 2. );
    breitWigner::resonance::Type
      base{ resonanceEnergy,
            neutronWidth,
            captureWidth,
            fissionWidth,
            weightedCompetitiveWidth,
            inversePenetrationFactor,
            shiftFactor,
            statisticalFactor };

    const auto& resonance =
      static_cast< const breitWigner::singleLevel::Resonance& >( base );
    
    auto resonanceXS = evaluate( resonance, resonanceEnergy ).data;
    auto resonanceElastic = std::get<0>(resonanceXS);
    auto resonanceCapture = std::get<1>(resonanceXS);
    auto resonanceFission = std::get<2>(resonanceXS);
    auto resonanceTotal = resonanceElastic
                          + resonanceCapture
                          + resonanceFission;

    GIVEN("an evaluation at the right full width half max energy"){
      auto energy = resonanceEnergy + 0.5 * totalWidth;
      
      auto fwhmXS = evaluate( resonance, energy ).data;
      auto fwhmElastic = std::get<0>( fwhmXS );
      auto fwhmCapture = std::get<1>( fwhmXS );
      auto fwhmFission = std::get<2>( fwhmXS );
      auto fwhmTotal = fwhmElastic
                       + fwhmCapture
                       + fwhmFission;

      THEN("the shape will be correct"){
        REQUIRE( static_cast< double >( 0.5 * resonanceTotal ) == 
                 Approx( static_cast< double >( fwhmTotal ) ) );
      }
      
      THEN("the magnitude will match an alternative formulation"){
        auto ps = penetrationShift( l, k( energy ) * a( energy ) );
      
        auto primedResonanceEnergy =
          resonance.energy
          + 0.5 * resonance.neutronWidth
            * resonance.inversePenetrationFactor
            * ( resonance.shiftFactor - ps[1] );

        auto neutronWidth =
          resonance.neutronWidth
          * ps[0]
          * resonance.inversePenetrationFactor;      
      
        // for some reason, gcc cannot deduce that reference and fwhmCapture
        // have the same unit
        // auto reference =
        auto reference = fwhmCapture;
        reference =
          resonance.statisticalFactor
          * neutronWidth
          * resonance.captureWidth
          / ( pow( resonance.captureWidth
                   + resonance.fissionWidth
                   + neutronWidth, Ratio<2> )
              + 4. * pow( energy - primedResonanceEnergy, Ratio<2> ) );

        REQUIRE( static_cast< double >( reference ) == 
                 Approx( static_cast< double >( fwhmCapture ) ) );
      }
    }
    
    GIVEN("an evaluation at the left full width half max energy") {
      auto energy = resonanceEnergy - 0.5 * totalWidth;
      auto fwhmXS = evaluate( resonance, energy ).data;
      auto fwhmElastic = std::get<0>(fwhmXS);
      auto fwhmCapture = std::get<1>(fwhmXS);
      auto fwhmFission = std::get<2>(fwhmXS);
      auto fwhmTotal = fwhmElastic
                       + fwhmCapture
                       + fwhmFission;

      THEN("the shape will be correct"){
        REQUIRE( static_cast< double >( 0.5 * resonanceTotal ) == 
                 Approx( static_cast< double >( fwhmTotal ) ) );
      }
      
      THEN("the magnitude will match an alternative formulation"){
        auto ps = penetrationShift( l, k( energy ) * a( energy ) );
      
        auto primedResonanceEnergy =
          resonance.energy
          + 0.5 * resonance.neutronWidth
            * resonance.inversePenetrationFactor
            * ( resonance.shiftFactor - ps[1] );

        auto neutronWidth =
          resonance.neutronWidth
          * ps[0]
          * resonance.inversePenetrationFactor;      
      
        // for some reason, gcc cannot deduce that reference and fwhmCapture
        // have the same unit
        // auto reference =
        auto reference = fwhmCapture;
        reference =
          resonance.statisticalFactor
          * neutronWidth
          * resonance.captureWidth
          / ( pow( resonance.captureWidth
                   + resonance.fissionWidth
                   + neutronWidth, Ratio<2> )
              + 4. * pow( energy - primedResonanceEnergy, Ratio<2> ) );

        REQUIRE( static_cast< double >( reference ) == 
                 Approx( static_cast< double >( fwhmCapture ) ) );
      }
    }    
  }
}
