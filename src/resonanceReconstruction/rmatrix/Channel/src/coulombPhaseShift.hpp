/**
 *  @brief Return the coulomb phase shift for this channel as a function of
 *         energy
 *
 *  @param[in] energy   the energy at which the coulomb phase shift is needed
 */
double coulombPhaseShift( const Energy& energy ) const {
  auto function = [&] ( auto type ) {
    const double eta = this->particlePair().etaParameter( energy ).value;
    const unsigned int l = this->quantumNumbers().orbitalAngularMomentum();
    return calculateCoulombPhaseShift< decltype( type ) >( l, eta );
  };
  return std::visit( function, this->type_ );
}

