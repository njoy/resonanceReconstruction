/**
 *  @brief Generate the R-matrix contribution for the resonance
 *
 *  The R-matrix element R_{ij} for channel i to i is given by:
 *      R_{ij} = \frac{ \gamma_i \gamma_j }{ E_r - E - i \Gamma/2 }
 *  in which \gamma_i and \gamma_j are the reduced widths of channel i and j,
 *  E_r is the resonance energy, E is the energy at which the resonance needs
 *  to be evaluated and \Gamma is the eliminated capture width.
 *
 *  Because the eliminated channel is a capture channel, the value of \Gamma
 *  is energy dependent and does not require a penetrability value.
 *
 *  @param[in] energy   the energy at which the resonance must be evaluated
 */
auto rmatrix( const Energy& energy ) const {

  const auto terminator = [&] { 
    const auto deltaEnergy = this->energy() - energy;
    const auto eliminatedWidth = this->eliminatedWidth();
    const auto halfCaptureWidth = eliminatedWidth * eliminatedWidth;
    return 1. / std::complex< double >(  deltaEnergy / electronVolt,
                                        -halfCaptureWidth / electronVolt ); 
  }();

  const auto widths = this->widths();
  const int size = widths.size();
  return ranges::view::cartesian_product( widths, widths )
           | ranges::view::transform(
                 [] ( auto&& pair ) -> Energy
                    { return std::get< 0 >( pair ) *
                             std::get< 1 >( pair ); } )
           | ranges::view::transform(
                 [] ( auto&& square ) -> double
                    { return square / electronVolt; } )
           | ranges::view::transform(
                 [terminator] ( auto&& element ) -> std::complex< double >
                              { return element * terminator; } )
           | ranges::view::chunk( size );
}

