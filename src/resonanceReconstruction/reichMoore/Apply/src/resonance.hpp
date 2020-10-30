template< typename ChannelRatio,
          typename StatisticalFactor,
          typename PenetrationShift,
          typename Range >
static auto
resonance( const endf::ReichMoore::LValue::Resonance< Range >& resonance,
           ChannelRatio&& rho,
           StatisticalFactor&& g,
           PenetrationShift&& penetrationShift ){
  const auto energy = resonance.ER() * electronVolts;
  const auto neutronWidth = resonance.GN() * electronVolts;
  const auto captureWidth = resonance.GG() * electronVolts;
  const auto fissionWidthA = resonance.GFA() * electronVolts;
  const auto fissionWidthB = resonance.GFB() * electronVolts;

  const double penetrationFactor = penetrationShift( rho( std::abs(energy) ) )[0];

  const double flaggedStatisticalFactor =
    std::copysign( g( std::abs( resonance.AJ() ) ), resonance.AJ() );

  return Resonance( energy,
                    neutronWidth,
                    captureWidth,
                    fissionWidthA,
                    fissionWidthB,
                    penetrationFactor,
                    flaggedStatisticalFactor );
}

template< typename Range, typename Factory >
static auto resonances( Range&& rs, Factory&& factory ){
  return rs | ranges::view::transform( factory ) | ranges::to_vector;
}
