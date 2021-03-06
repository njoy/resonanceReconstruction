/**
 *  @brief Return the phase shift for this channel as a function of energy
 *
 *  @param[in] energy   the energy at which the phase shift is needed
 */
double phaseShift( const Energy& energy ) const {

  const double eta = this->sommerfeldParameter( energy );
  const double ratio = this->waveNumber( energy ) *
                       this->radii().phaseShiftRadius( energy );
  const unsigned int l = this->quantumNumbers().orbitalAngularMomentum();
  return calculatePhaseShift< ChannelType >( l, ratio, eta );
}
