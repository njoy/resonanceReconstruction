/**
 *  @brief Calculate the ( I - RL )^-1 R matrix
 *
 *  @param[in] energy            the energy value
 *  @param[in] table             the resonance table
 *  @param[in] penetrabilities   the channel penetrabilities
 *  @param[in] channels          the channels
 *  @param[in] belowThreshold    which channels are below the threshold
 *
 *  @return The resulting ( I - RL )^-1 R matrix
 */
template < typename Penetrabilities, typename Channels, typename BelowThreshold >
const Matrix< std::complex< double > >&
operator()( const Energy& energy,
            const ResonanceTable& table,
            const Penetrabilities& penetrabilities,
            const Channels& channels,
            const BelowThreshold& belowThreshold ) {

  //! @todo handle threshold reactions?

  // calculate the L matrix
  this->lmatrix_( energy, penetrabilities, channels );

  // calculate the level matrix
  auto energies = table.energies();
  const unsigned int numberResonances = table.numberResonances();
  const unsigned int numberChannels = table.numberChannels();
  this->amatrix_ = Matrix< std::complex< double > >::Zero( numberResonances,
                                                           numberResonances );
  for ( unsigned int lambda = 0; lambda < numberResonances; ++lambda ) {
    this->amatrix_( lambda, lambda ) = energies[lambda] - energy;
    for ( unsigned int mu = 0; mu < numberResonances; ++mu ) {
      for ( unsigned int c = 0; c < numberChannels; ++c ) {
        this->amatrix_( lambda, mu ) -= this->gmatrix_( c, lambda ) *
                                        this->lmatrix_.diagonal()[c] *
                                        this->gmatrix_( c, mu ) /
                                        electronVolt;
      }
    }
  }
  this->amatrix_ = this->amatrix_.inverse();

  // calculate and return R_L = ( 1 - RL )^-1 R = G A G^T
  this->rlmatrix_ = this->gmatrix_ * this->amatrix_ *
                    this->gmatrix_.transpose();
  return this->rlmatrix_;
}
