void precompute( const Energy& energy ) {

  // data we need: P(E), P(E-Q), S(E)
  const auto channel = this->incidentChannel();
  const auto qx = this->QX();
  const auto p = channel.penetrability( energy );
  const auto q = qx.value != 0. ? channel.penetrability( energy - qx )
                                : p;
  const auto s = channel.shiftFactor( energy );

  // useful lambdas
  auto elastic = [&] ( const auto& resonance ) { return resonance.elastic( p ); };
  auto capture = [] ( const auto& resonance ) { return resonance.capture(); };
  auto fission = [] ( const auto& resonance ) { return resonance.fission(); };
  auto total = [&] ( const auto& resonance ) { return resonance.total( p, q ); };
  auto delta = [&] ( const auto& resonance )
                   { return energy - resonance.energyPrime( s ); };
  auto denominator = [] ( const auto& delta, const auto& total )
                        { return delta * delta + 0.25 * total * total; };

  // precompute values
  const auto resonances = this->resonanceTable().resonances();
  this->elastic_ = resonances | ranges::view::transform( elastic );
  this->capture_ = resonances | ranges::view::transform( capture );
  this->fission_ = resonances | ranges::view::transform( fission );
  this->total_ = resonances | ranges::view::transform( total );
  this->delta_ = resonances | ranges::view::transform( delta );
  this->denominator_ = ranges::view::zip_with( denominator,
                                               this->delta_, this->total_ );
}
