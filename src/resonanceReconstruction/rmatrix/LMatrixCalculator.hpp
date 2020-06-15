template < typename BoundaryOption >
class LMatrixCalculator;

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
    lmatrix_( numberChannels ) {};

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

    auto real = [&] ( const auto& channel )
                    { return channel.shiftFactor( energy ) -
                             channel.boundaryCondition(); };
    auto complex = [] ( double real, double imaginary )
                      { return std::complex< double >( real, imaginary ); };

    auto diagonal = ranges::view::zip_with(
                      complex,
                      channels | ranges::view::transform(
                                   [&] ( const auto& channel )
                                       { return std::visit( real, channel ); } ),
                      penetrabilities );

    this->lmatrix_.setZero();
    for ( unsigned int i = 0; i < diagonal.size(); ++i ) {

      this->lmatrix_.diagonal()[i] = diagonal[i];
    }
    return this->lmatrix_;
  }
};

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
    lmatrix_( numberChannels ) {};

  /**
  *  @brief Diagonal L matrix calculation for the ShiftFactor boundary condition
  *
  *  @param[in] energy            the energy value (unused)
  *  @param[in] penetrabilities   the channel penetrabilities
  *  @param[in] channels          the channels (unused)
  *
  *  @return The resulting L = S - B + iP matrix
   */
  template < typename Penetrabilities, typename Channels >
  const DiagonalMatrix< std::complex< double > >&
  operator()( const Energy&, const Penetrabilities& penetrabilities,
              const Channels& ) {

    this->lmatrix_.setZero();
    for ( unsigned int i = 0; i < penetrabilities.size(); ++i ) {

      this->lmatrix_.diagonal()[i] = std::complex< double >( 0.0, penetrabilities[i] );
    }
    return this->lmatrix_;
  }
};
