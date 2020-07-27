/**
 *  @brief Evaluate the cross sections at the given energy
 *
 *  @param[in] energy       the incident energy
 *  @param[in,out] result   a map containing the accumulated cross sections
 */
void evaluate( const Energy& energy,
               std::map< ReactionID, CrossSection >& result ) {

// ----- DEBUG -----
std::cout << "-- CompoundSystem is working on this energy: " << energy << " --" << std::endl;
// ----- DEBUG -----

  // accumulate over each spin group
  ranges::for_each( this->groups_,
                    [&] ( auto& group )
                        { group.evaluate( energy, result ); } );

  // calculate potential scattering
  const auto channel = this->groups_.front().incidentChannel();
  const auto incident = channel.particlePair().pairID();
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

// ----- DEBUG -----
bool fission = ( result.find( ReactionID( incident.symbol() + "->fission" ) ) != result.end() );
std::cout << "potential " << factor * value << std::endl;
  result[ ReactionID( incident, incident ) ] += factor * value;
std::cout << "final " << result[ ReactionID( incident, incident ) ] << " "
                      << result[ ReactionID( incident.symbol() + "->capture" ) ] << " "
                      << ( fission == true ? result[ ReactionID( incident.symbol() + "->fission" ) ] : 0.*barns ) << std::endl;
std::cout << "-- CompoundSystem is done for this energy --" << std::endl;
// ----- DEBUG -----
}