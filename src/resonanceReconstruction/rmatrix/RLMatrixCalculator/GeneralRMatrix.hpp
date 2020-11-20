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
   *  @brief Return the L = S - B + iP matrix
   */
  const DiagonalMatrix< std::complex< double > >&
  lmatrix() const { return this->lmatrix_.matrix(); }

  /**
   *  @brief Return the G matrix
   */
  const Matrix< std::complex< double > >&
  gmatrix() const { return this->gmatrix_; }

  /**
   *  @brief Return the A matrix
   */
  const Matrix< std::complex< double > >&
  amatrix() const { return this->amatrix_; }

  /**
   *  @brief Return the ( I - RL )^-1 R matrix
   */
  const Matrix< std::complex< double > >&
  rlmatrix() const { return this->rlmatrix_; }

  #include "resonanceReconstruction/rmatrix/RLMatrixCalculator/GeneralRMatrix/src/call.hpp"
};
