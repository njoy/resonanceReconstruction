/**
 *  @brief Return the statistical spin factor
 *
 *  The statistical spin factor g of a channel is defined as follows:
 *     g = ( 2 * J + 1 ) / ( 2 * ia + 1 ) / ( 2 * ib + 1 )
 *  in which J is the total angular momentum of the channel and ia and ib 
 *  are the spins of the particles in the particle pair.
 *
 *  For a particle pair involving a neutron (for which the particle spin is 
 *  0.5), this reduces to:
 *    g = ( 2 * J + 1 ) / ( 2 * I + 1 ) / 2
 *  in which I is the target nucleus spin value.
 */
double statisticalSpinFactor() const {

  auto J = this->quantumNumbers().totalAngularMomentum();
  auto ia = this->particlePair().particle().spin();
  auto ib = this->particlePair().residual().spin();
  return  ( 2. * J + 1. ) / ( 2. * ia + 1. ) / ( 2. * ib + 1. );
}

