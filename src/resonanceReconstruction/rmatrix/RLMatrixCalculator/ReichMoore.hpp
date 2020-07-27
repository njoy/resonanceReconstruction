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
   *  @brief Return the ( I - RL )^-1 R matrix
   */
  const Matrix< std::complex< double > >&
  matrix() const { return this->lmatrix_; }

  #include "resonanceReconstruction/rmatrix/RLMatrixCalculator/ReichMoore/src/call.hpp"
};
