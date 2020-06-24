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

  this->lmatrix_.setZero();
  for ( unsigned int i = 0; i < penetrabilities.size(); ++i ) {

    this->lmatrix_.diagonal()[i] =
        std::complex< double >( 0.0, penetrabilities[i] );
  }
  return this->lmatrix_;
}
