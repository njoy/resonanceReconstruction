/**
 *  @brief Return the value of the Mandelstam variable s parameter for the
 *         particle pair at a given energy
 *
 *  The Mandelstam variable s is defined as:
 *     s = ( ( ma + mb ) c^2 )^2 + 2 mb c^2 energy
 *  in which ma and mb are the atomic masses of the first and second particles
 *  in the particle pair and c is the light speed (converting the mass into
 *  energy).
 *
 *  @param energy   the energy at which the Mandelstam variable s needs to be
 *                  evaluated
 */
EnergySquared mandelstam( const Energy& energy ) const {

  const auto particle = this->particlePair().particle().mass() * c * c;
  const auto residual = this->particlePair().residual().mass() * c * c;
  const auto pair = particle + residual;
  return pair * pair + 2. * residual * energy;
}
