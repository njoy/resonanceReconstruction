/**
 *  @class
 *  @brief A functor to calculate the L matrix using the ShiftFactor boundary
 *         condition
 *
 *  The L matrix is a diagonal matrix defined as S - B + iP
 *
 *  The R-matrix analysis code SAMMY uses a different approach to the boundary
 *  condition B. SAMMY effectively eliminates the real part of the L matrix
 *  diagonal by setting the boundary condition B to be equal to the shift
 *  factor S, making the bondary condition potentially an energy dependent
 *  quantity.
 *
 *  Using the SAMMY boundary condition will effectively ignore the value of the
 *  boundary condition of each channel.
 */
template <>
class LMatrixCalculator< ShiftFactor > {

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

  #include "resonanceReconstruction/rmatrix/LMatrixCalculator/ShiftFactor/src/call.hpp"
};
