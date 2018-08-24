/**
 *  @brief Return the shift factor for this channel as a function of energy
 *
 *  @param[in] energy   the energy at which the shift factor is needed
 */
double shiftFactor( const Energy& energy ) const {

  const double eta = this->particlePair().etaParameter( energy ).value;
  const double ratio = this->particlePair().waveNumber( energy ) *
                       this->radii().shiftFactorRadius( energy );
  const unsigned int l = this->quantumNumbers().orbitalAngularMomentum();
  return calculateShiftFactor< ChannelType >( l, ratio, eta );
}

