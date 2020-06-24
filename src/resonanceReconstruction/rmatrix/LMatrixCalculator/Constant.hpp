/**
 *  @class
 *  @brief A functor to calculate the diagonal L matrix using the Constant
 *         boundary condition
 *
 *  The L matrix is a diagonal matrix defined as S - B + iP
 */
template <>
class LMatrixCalculator< Constant > {

  /* fields */
  DiagonalMatrix< std::complex< double > > lmatrix_;

public:

  /**
   *  @brief Constructor
   *
   *  @param[in] numberChannels   the number of channels
   */
  LMatrixCalculator( unsigned int numberChannels ) :
    lmatrix_( numberChannels ) {

    this->lmatrix_.setZero();
  };

  /**
   *  @brief Return the diagonal L matrix
   */
  const DiagonalMatrix< std::complex< double > >&
  matrix() const { return this->lmatrix_; }

  #include "resonanceReconstruction/rmatrix/LMatrixCalculator/Constant/src/call.hpp"
};
