/**
 *  @brief Evaluate the cross sections at the given energy
 *
 *  @param[in] energy       the incident energy
 *  @param[in,out] result   a map containing the accumulated cross sections
 */
void evaluate( const Energy& energy,
               Map< ReactionID, CrossSection >& result ) {

  // accumulate over each spin group
  ranges::for_each( this->groups_,
                    [&] ( auto& group )
                        { group.evaluate( energy, result ); } );

  // calculate potential scattering
  const auto channel = this->groups_.front().incidentChannel();
  const auto incident = channel.particlePair().particle().particleID();
  const auto target = channel.particlePair().residual().particleID();
  const auto waveNumber = channel.waveNumber( energy );
  const auto ratio = waveNumber * channel.radii().phaseShiftRadius( energy );

  // the 4 * pi / k2 factor
  const CrossSection factor =  4. * pi / ( waveNumber * waveNumber );

  // accumulate potential scattering
  double value = 0;
  for ( unsigned int l = 0; l <= this->lmax_; ++l ) {

    const double phi = calculatePhaseShift< Neutron >( l, ratio, 0. );
    const double sinphi = std::sin( phi );
    const double sin2phi = sinphi * sinphi;
    value += ( 2. * l + 1. ) * sin2phi;
  }

  result[ ReactionID( incident, target, elementary::ReactionType( "elastic" ) ) ] += factor * value;
}
