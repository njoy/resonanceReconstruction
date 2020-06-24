/**
 *  @brief Diagonal L matrix calculation for the Constant boundary condition
 *
 *  @param[in] energy            the energy value
 *  @param[in] penetrabilities   the channel penetrabilities
 *  @param[in] channels          the channels
 *
 *  @return The resulting L = S - B + iP matrix
 */
template < typename Penetrabilities, typename Channels >
const DiagonalMatrix< std::complex< double > >&
operator()( const Energy& energy,
            const Penetrabilities& penetrabilities,
            const Channels& channels ) {

  auto shiftMinusBoundary = [&] ( const auto& channel )
                                { return channel.shiftFactor( energy ) -
                                         channel.boundaryCondition(); };
  auto toComplex = [] ( double real, double imaginary )
                      { return std::complex< double >( real, imaginary ); };

  auto diagonal = ranges::view::zip_with(
                    toComplex,
                    channels | ranges::view::transform(
                                 [&] ( const auto& channel )
                                     { return std::visit( shiftMinusBoundary,
                                                          channel ); } ),
                    penetrabilities );

  this->lmatrix_.setZero();
  for ( unsigned int i = 0; i < diagonal.size(); ++i ) {

    this->lmatrix_.diagonal()[i] = diagonal[i];
  }
  return this->lmatrix_;
}
