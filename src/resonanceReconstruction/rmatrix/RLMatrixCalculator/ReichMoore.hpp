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
    rlmatrix_( table.numberChannels(), table.numberChannels() ) {

    this->rmatrix_.setZero();
    this->rlmatrix_.setZero();
  };

  /**
   *  @brief Return the L = S - B + iP matrix
   */
  const DiagonalMatrix< std::complex< double > >&
  lmatrix() const { return this->lmatrix_.matrix(); }

  /**
   *  @brief Return the R matrix
   */
  const Matrix< std::complex< double > >&
  rmatrix() const { return this->rmatrix_; }

  /**
   *  @brief Return the ( I - RL )^-1 R matrix
   */
  const Matrix< std::complex< double > >&
  rlmatrix() const { return this->rlmatrix_; }

  #include "resonanceReconstruction/rmatrix/RLMatrixCalculator/ReichMoore/src/call.hpp"
};
