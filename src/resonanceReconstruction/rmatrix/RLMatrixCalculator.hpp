template < typename Formalism, typename BoundaryOption >
class RLMatrixCalculator;

/**
 *  @class
 *  @brief A functor to calculate the ( I - RL )^-1 R matrix using the
 *         Reich-Moore approximation
 */
template < typename BoundaryOption >
class RLMatrixCalculator< ReichMoore, BoundaryOption > {

  /* fields */
  LMatrixCalculator< BoundaryOption > lmatrix_;
  Matrix< std::complex< double > > rmatrix_;
  Matrix< std::complex< double > > rlmatrix_;

public:

  /**
   *  @brief Constructor
   *
   *  @param[in] table   the resonance table
   */
  RLMatrixCalculator( const ResonanceTable& table ) :
    lmatrix_( table.numberChannels() ),
    rmatrix_( table.numberChannels(), table.numberChannels() ),
    rlmatrix_( table.numberChannels(), table.numberChannels() ) {};

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

    // range with the R-matrices for each resonance
    auto rmatrices = table.resonances()
                       | ranges::view::transform(
                           [&] ( const auto& resonance )
                               { return resonance.rmatrix( energy ); } );

    // accumulate the rmatrix
    const unsigned int size = table.numberChannels();
    this->rmatrix_ = Matrix< std::complex< double > >::Zero( size, size );
    for ( const auto& rmatrix : rmatrices ) {
      for ( unsigned int c = 0; c < size; ++c ) {
        for ( unsigned int cprime = 0; cprime < size; ++cprime ) {
          this->rmatrix_( c, cprime ) += rmatrix[c][cprime];
        }
      }
    }

    // zero out threshold reactions
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
};
