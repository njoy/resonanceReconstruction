static auto competitiveWidth( const Quantity< ElectronVolts > offset ){
  return
    [offset]( auto&& resonance, auto&& rho, auto&& penetrationShift ) {

    const auto energy = resonance.ER() * electronVolts + offset;
    const auto competitiveWidth = resonance.GX() * electronVolts;

    return competitiveWidth / penetrationShift( rho( energy ) )[0];
  };
}
