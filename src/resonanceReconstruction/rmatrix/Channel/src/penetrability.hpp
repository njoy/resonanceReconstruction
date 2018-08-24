/**
 *  @brief Return the penetrability for this channel as a function of energy
 *
 *  @param[in] energy   the energy at which the penetrability is needed
 */
double penetrability( const Energy& energy ) const {

  const double eta = this->particlePair().etaParameter( energy ).value;
  const double ratio = this->particlePair().waveNumber( energy ) *
                       this->radii().penetrabilityRadius( energy );
  const unsigned int l = this->quantumNumbers().orbitalAngularMomentum();
  return calculatePenetrability< ChannelType >( l, ratio, eta );
}

