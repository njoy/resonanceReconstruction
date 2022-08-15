void precompute( const Energy& energy ) {

  // data we need: P(E), P(E-Q), S(E)
  const auto channel = this->incidentChannel();
  const auto qx = this->QX();
  const auto p = channel.penetrability( energy );
  const auto q = qx.value != 0. ? channel.penetrability( energy - qx )
                                : p;
  const auto s = channel.shiftFactor( energy );

  // useful lambdas
  auto elastic = [&] ( const auto& resonance ) -> decltype(auto) {

    return resonance.elastic( p );
  };
  auto capture = [] ( const auto& resonance ) -> decltype(auto) {

    return resonance.capture();
  };
  auto fission = [] ( const auto& resonance ) -> decltype(auto) {

    return resonance.fission();
  };
  auto total = [&] ( const auto& resonance ) -> decltype(auto) {

    return resonance.total( p, q );
  };
  auto delta = [&] ( const auto& resonance ) -> Energy {

    return energy - resonance.energyPrime( s );
  };
  auto denominator = [] ( const auto& delta, const auto& total ) -> EnergySquared { 

    return delta * delta + 0.25 * total * total;
  };

  // precompute values
  const auto resonances = this->resonanceTable().resonances();
  this->elastic_ = ranges::to< std::vector< Width > >(
                       resonances | ranges::cpp20::views::transform( elastic ) );
  this->capture_ = ranges::to< std::vector< Width > >(
                       resonances | ranges::cpp20::views::transform( capture ) );
  this->fission_ = ranges::to< std::vector< Width > >(
                       resonances | ranges::cpp20::views::transform( fission ) );
  this->total_ = ranges::to< std::vector< Width > >(
                     resonances | ranges::cpp20::views::transform( total ) );
  this->delta_ = ranges::to< std::vector< Energy > >(
                     resonances | ranges::cpp20::views::transform( delta ) );
  this->denominator_ = ranges::to< std::vector< EnergySquared > >(
                           ranges::views::zip_with( denominator,
                                                    this->delta_,
                                                    this->total_ ) );
}
