/**
 *  @brief Return the relativistic wave number of the particle pair at a given
 *         energy
 *
 *  The relativistic wave number k is an energy dependent quantity defined as
 *  follows:
 *     k^2 = ( s - ( ( ma + mb ) c^2 )^2 ) ( s - ( ( ma - mb ) c^2 )^2 ) / 4 / s
 *  in which s is the Mandelstam variable s, ma and mb are the atomic masses
 *  of the first and second particles in the particle pair and c is the light
 *  speed (converting the mass into energy).
 *
 *  @param energy   the energy at which the relativistic wave number needs
 *                  to be evaluated
 */
WaveNumber waveNumber( const Energy& energy ) const {

  const auto s = this->mandelstam( energy );
  const auto particle = this->particle().mass() * c * c;
  const auto residual = this->residual().mass() * c * c;
  const auto pair = particle + residual;
  const auto delta = particle - residual;
  return 0.5 * sqrt( ( s - pair * pair ) * ( s - delta * delta ) / s );
}
