template< typename ChannelRatio,
          typename StatisticalFactor,
          typename PenetrationShift >
static auto resonance( const ENDF::resolved::SLBW::Lvalue::Resonance& resonance,
                       ChannelRatio&& rho,
                       StatisticalFactor&& g,
                       PenetrationShift&& penetrationShift ){
  const auto energy = resonance.ER() * electronVolts;
  const auto neutronWidth = resonance.GN() * electronVolts;
  const auto captureWidth = resonance.GG() * electronVolts;
  const auto fissionWidth = resonance.GF() * electronVolts;

  const auto competitiveWidth = 0.0 * electronVolts;

  const auto ps = penetrationShift( rho( energy ) );
  const auto& penetrationFactor = ps[0];
  const auto& shiftFactor = ps[1];

  const auto statisticalFactor = g( std::abs( resonance.AJ() ) );

  return resonance::Type( energy,
                          neutronWidth,
                          captureWidth,
                          fissionWidth,
                          competitiveWidth,
                          penetrationFactor,
                          shiftFactor,
                          statisticalFactor );
}

template< typename ChannelRatio,
          typename StatisticalFactor,
          typename PenetrationShift,
          typename CompetitiveWidth >
static auto resonance( const ENDF::resolved::SLBW::Lvalue::Resonance& resonance,
                       ChannelRatio&& rho,
                       StatisticalFactor&& g,
                       PenetrationShift&& penetrationShift,
                       CompetitiveWidth&& GX ){
  const auto energy = resonance.ER() * electronVolts;
  const auto neutronWidth = resonance.GN() * electronVolts;
  const auto captureWidth = resonance.GG() * electronVolts;
  const auto fissionWidth = resonance.GF() * electronVolts;

  const auto competitiveWidth = GX( resonance, rho, penetrationShift );

  const auto ps = penetrationShift( rho( std::abs(energy) ) );
  const auto& penetrationFactor = ps[0];
  const auto& shiftFactor = ps[1];

  const auto statisticalFactor = g( std::abs( resonance.AJ() ) );

  return resonance::Type( energy,
                          neutronWidth,
                          captureWidth,
                          fissionWidth,
                          competitiveWidth,
                          penetrationFactor,
                          shiftFactor,
                          statisticalFactor );
}

template< typename Range, typename Factory >
static auto resonances( Range&& rs, Factory&& factory ){
  return rs | ranges::view::transform( factory ) | ranges::to_vector;
}
