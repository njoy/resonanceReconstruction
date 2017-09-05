struct ZeroWidth {
  template< typename Resonance >
  constexpr Quantity< ElectronVolts > operator()( const Resonance& ){
    return 0.0 * electronVolts;
  }
};
