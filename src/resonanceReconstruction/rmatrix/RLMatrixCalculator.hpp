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

/**
 *  @class
 *  @brief A functor to calculate the ( I - RL )^-1 R matrix in general
 *         R-matrix
 */
template < typename BoundaryOption >
class RLMatrixCalculator< GeneralRMatrix, BoundaryOption > {

  /* fields */
  LMatrixCalculator< BoundaryOption > lmatrix_;
  Matrix< std::complex< double > > gmatrix_;
  Matrix< std::complex< double > > amatrix_;
  Matrix< std::complex< double > > rlmatrix_;

  /* fields */
  #include "resonanceReconstruction/rmatrix/RLMatrixCalculator/src/makeGMatrix.hpp"

public:

  /**
   *  @brief Constructor
   *
   *  @param[in] table   the resonance table
   */
  RLMatrixCalculator( const ResonanceTable& table ) :
    lmatrix_( table.numberChannels() ),
    gmatrix_( makeGMatrix( table ) ),
    amatrix_( table.numberResonances(), table.numberResonances() ),
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
};
