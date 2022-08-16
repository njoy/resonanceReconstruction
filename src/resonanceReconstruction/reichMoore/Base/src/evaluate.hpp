auto evaluate( const Quantity<ElectronVolts> energy,
               const Quantity<RootBarn> channelRadius,
               const Quantity<RootBarn> scatteringRadius ) const {

  const Quantity< InvRootBarn > waveNumber = this->neutronWaveNumber( energy );
  const double channelRatio = waveNumber * channelRadius;
  const double scatteringRatio = waveNumber * scatteringRadius;

  auto dispatch = [&]( auto tag ){
    const auto crossSections =
      this->lvalues
      | ranges::views::transform
        ( Inspector< decltype( tag ) >
          { energy, waveNumber, channelRatio, scatteringRatio } );

    return ranges::accumulate( crossSections, pack( 0.0, 0.0, 0.0 ) ).data;
  };

  const auto reduction = std::visit( dispatch, this->tag );

  const auto scaling = pi / ( waveNumber * waveNumber );
  const auto elastic = std::get<0>( reduction ) * scaling;
  const auto capture = std::get<1>( reduction ) * scaling;
  const auto fission = std::get<2>( reduction ) * scaling;
  return CrossSection( elastic, capture, fission );
}
