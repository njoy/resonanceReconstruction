/**
 *  @brief Calculate the ( I - RL )^-1 R matrix
 *
 *  @param[in] energy            the energy value
 *  @param[in] table             the resonance table
 *  @param[in] penetrabilities   the channel penetrabilities
 *  @param[in] channels          the channels
 *
 *  @return The resulting ( I - RL )^-1 R matrix
 */
template < typename Penetrabilities, typename Channels >
const Matrix< std::complex< double > >&
operator()( const Energy& energy,
            const ResonanceTable& table,
            const Penetrabilities& penetrabilities,
            const Channels& channels ) {

  // accumulate the rmatrix
  const unsigned int size = table.numberChannels();
  this->rmatrix_ = Matrix< std::complex< double > >::Zero( size, size );
  for ( const auto& resonance : table.resonances() ) {

    decltype(auto) rmatrix = resonance.rmatrix( energy );
    for ( unsigned int c = 0; c < size; ++c ) {

      for ( unsigned int cprime = 0; cprime < size; ++cprime ) {

        this->rmatrix_( c, cprime ) += rmatrix[c][cprime];
      }
    }
  }

  // zero out threshold reactions
  auto belowThreshold = channels
         | ranges::views::transform(
              [&] ( const auto& channel )
                  { return std::visit(
                               [&] ( const auto& channel )
                                   { return channel.belowThreshold( energy ); },
                               channel ); } );
  for ( unsigned int c = 0; c < size; ++c ) {

    if ( belowThreshold[c] ) {

      this->rmatrix_.row(c).setZero();
      this->rmatrix_.col(c).setZero();
    }
  }

  // calculate and return R_L = ( 1 - RL )^-1 R
  this->rlmatrix_ = Matrix< double >::Identity( size, size );
  this->rlmatrix_ -= this->rmatrix_ *
                     this->lmatrix_( energy, penetrabilities, channels );
  this->rlmatrix_ = this->rlmatrix_.inverse();
  this->rlmatrix_ *= this->rmatrix_;
  return this->rlmatrix_;
}
