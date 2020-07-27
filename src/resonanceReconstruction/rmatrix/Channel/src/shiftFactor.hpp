/**
 *  @brief Return the shift factor for this channel as a function of energy
 *
 *  @param[in] energy   the energy at which the shift factor is needed
 */
double shiftFactor( const Energy& energy ) const {

  const double eta = this->sommerfeldParameter( energy );
  const double ratio = this->waveNumber( energy ) *
                       this->radii().shiftFactorRadius( energy );
  const unsigned int l = this->quantumNumbers().orbitalAngularMomentum();
  return calculateShiftFactor< ChannelType >( l, ratio, eta );
}
