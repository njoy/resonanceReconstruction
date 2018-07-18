/**
 *  @brief Return the reduced mass of the particle pair
 *
 *  The reduced mass mu of the two particles is defined as follows:
 *     mu = ma * mb / ( ma + mb )
 *  in which ma and mb are the atomic mass values of the particles in the
 *  particle pair.
 */
AtomicMass reducedMass() const {

  auto ma = this->particle().mass();
  auto mb = this->residual().mass();
  return ma * mb / ( ma + mb );
}

