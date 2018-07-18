/**
 *  @brief Return the penetrability for this channel as a function of energy
 *
 *  @param[in] energy   the energy at which the penetrability is needed
 */
double penetrability( const Energy& energy ) const {

  auto function = [&] ( auto type ) {
    const double ratio = this->particlePair().waveNumber( energy ) *
                         this->radii().penetrabilityRadius( energy );
    const unsigned int l = this->quantumNumbers().orbitalAngularMomentum();
    return calculatePenetrability< decltype( type ) >( l, ratio );
  };
  return std::visit( function, this->type_ );
}

