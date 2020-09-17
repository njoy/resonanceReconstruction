template < typename LValue >
std::vector< legacy::unresolved::SpinGroup >
makeLegacyUnresolvedSpinGroups(
    const LValue& endfLValue,
    const ParticlePair& pair, const ChannelRadii& radii,
    const std::vector< double >& energies ) {

  // useful numbers
  unsigned int l = endfLValue.orbitalMomentum();
  unsigned int njs = endfLValue.numberSpinValues();

  // some useful lambdas
  auto makeSpinGroup = [&] ( unsigned int l, const auto& in, const auto& radii,
                             const auto& endfJValue ) {

    return legacy::unresolved::SpinGroup(
               Channel< Neutron >( in, in, 0. * electronVolt,
                                   { l, 0.5, endfJValue.spin(),
                                     std::pow( -1, l ) > 0.
                                         ? static_cast< Parity >( +1 )
                                         : static_cast< Parity >( -1 ) },
                                   radii ),
               makeLegacyUnresolvedResonanceTable( endfJValue, energies ) );
  };

  // some ranges magic
  auto repeatL = ranges::view::repeat_n( l, njs );
  auto repeatPair = ranges::view::repeat_n( pair, njs );
  auto repeatRadii = ranges::view::repeat_n( radii, njs );
  auto jvalues = endfLValue.jValues();

  return ranges::view::zip_with(
           makeSpinGroup,
           repeatL,
           repeatPair,
           repeatRadii,
           jvalues );
}
