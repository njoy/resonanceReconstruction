static auto competitiveWidth( const Quantity< ElectronVolts > offset ){
  return
    [offset]( const ENDF::resolved::MLBW::LState::Resonance& resonance,
              auto&& rho, auto&& penetrationShift ){
    
    const auto energy = resonance.ER() * electronVolts + offset;
    const auto competitiveWidth = ( resonance.GT() 
                                    - resonance.GN() 
                                    - resonance.GG() 
                                    - resonance.GF() ) * electronVolts;

    return competitiveWidth / penetrationShift( rho( energy ) )[0];
  };
}
