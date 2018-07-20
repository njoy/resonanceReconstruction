/**
 *  @class
 *  @brief A spin group corresponding to a J,pi value
 */
class SpinGroup {

  /* fields */
  std::vector< Channel > channels_; // first channel is entrance
  ResonanceTable parameters_;

  auto penetrabilities( const Energy& energy ) const {
    return this->channels_
             | ranges::view::transform(
                 [&] ( const auto& channel )
                     { return channel.penetrability( energy ); } );
  }

  auto shiftFactors( const Energy& energy ) const {
    return this->channels_
             | ranges::view::transform(
                 [&] ( const auto& channel )
                     { return channel.shiftFactor( energy ); } );
  }

  auto phaseShifts( const Energy& energy ) const {
    return this->channels_
             | ranges::view::transform(
                 [&] ( const auto& channel )
                     { return channel.phaseShift( energy ); } );
  }

  auto coulombShifts( const Energy& energy ) const {
    return this->channels_
             | ranges::view::transform(
                 [&] ( const auto& channel )
                     { return channel.coulombPhaseShift( energy ); } );
  }

  auto boundaryConditions() const {
    return this->channels_
             | ranges::view::transform(
                 [&] ( const auto& channel )
                     { return channel.boundaryCondition(); } );
  }

public:

  /* constructor */
  SpinGroup( std::vector< Channel >&& channels, ResonanceTable&& table ) :
    channels_( std::move( channels ) ), parameters_( std::move( table ) ) {}

  auto resonanceTable() const { return this->parameters_; }
  auto channels() const { return ranges::view::all( this->channels_ ); }
  auto incidentChannel() const { return this->channels_.front(); }

  auto evaluate( const Energy& energy ) const {

    // penetrabilities, shift factors, phase shifts and Coulomb phase shifts
    // for each channel
    const auto penetrabilities = this->penetrabilities( energy );
    const auto shiftFactors = this->shiftFactors( energy );
    const auto phaseShifts = this->phaseShifts( energy );
    const auto coulombShifts = this->coulombShifts( energy );
    const auto boundaryConditions = this->boundaryConditions();

std::cout << "penetrabilities(0) " << penetrabilities[0] << std::endl;
std::cout << "shiftFactors(0) " << shiftFactors[0] << std::endl;
std::cout << "phaseShifts(0) " << phaseShifts[0] << std::endl;
std::cout << "coulombShifts(0) " << coulombShifts[0] << std::endl;
std::cout << "boundaryConditions(0) " << boundaryConditions[0] << std::endl;

    // the diagonal of the L matrix = S - B, iP
    const auto diagonalLMatrix =
      ranges::view::zip_with(
        [] ( const double a, const double b )
           { return std::complex< double >( a, b ); },
        ranges::view::zip_with( std::minus<double>(),
                                shiftFactors, boundaryConditions ),
        penetrabilities );

std::cout << "diagonalLMatrix(0) " << diagonalLMatrix[0] << std::endl;

    // the diagonal of the sqrt(P) matrix
    const auto diagonalSqrtPMatrix =
      penetrabilities
        | ranges::view::transform(
            [] ( const double penetrability ) -> double
               { return std::sqrt( penetrability ); } );

std::cout << "diagonalSqrtPMatrix(0) " << diagonalSqrtPMatrix[0] << std::endl;

    // the X matrix multiplier for each channel
    const auto xMultipliers =
      diagonalSqrtPMatrix
        | ranges::view::transform(
            [&] ( const double value ) -> double
                { return diagonalSqrtPMatrix.front() * value; } );

std::cout << "xMultipliers(0) " << xMultipliers[0] << std::endl;

    // the diagonal of the Omega matrix = exp( i ( w - phi ) )
    const auto diagonalOmegaMatrix =
      ranges::view::zip_with( std::minus< double >(),
                              coulombShifts, phaseShifts )
        | ranges::view::transform(
            [] ( const double delta ) -> std::complex< double >
               { return std::exp( std::complex< double >( 0.0, delta ) ); } );

std::cout << "diagonalOmegaMatrix(0) " << diagonalOmegaMatrix[0] << std::endl;

    // the U matrix multiplier for each channel
    const auto uMultipliers =
      diagonalOmegaMatrix
        | ranges::view::transform(
            [&] ( const std::complex< double > value ) -> std::complex< double >
                { return diagonalOmegaMatrix.front() * value; } );

std::cout << "uMultipliers(0) " << uMultipliers[0] << std::endl;

// BEGIN REALLY BAD - GET TESTING GOING

    //! @todo we only need the row of ( 1 - RL )^-1 corresponding to the 
    //!       incident channel and each column/line of R (R is a symmetric
    //!       matrix)

    // get the L matrix
    const unsigned int size = this->channels().size();
    Matrix< std::complex< double > > lMatrix =
        Matrix< std::complex< double > >::Zero( size, size );
    for ( unsigned int i = 0; i < size; ++i ) {
      lMatrix.diagonal()[i] = diagonalLMatrix[i];
    }

    // get the R-matrix and calculate ( 1 - RL )^-1 R
    Matrix< std::complex< double > > rMatrix =
        this->resonanceTable().rmatrix( energy );
    Matrix< std::complex< double > > xMatrix =
        ( Matrix< double >::Identity( size, size )
              - rMatrix * lMatrix ).inverse() * rMatrix;

    std::vector< std::complex<double> > row;
    for ( unsigned int i = 0; i < size; ++i ) row.push_back(xMatrix(0,i));

// END REALLY BAD - GET TESTING GOING

std::cout << "row(0) " << row[0] << std::endl;

    // the first row of the X matrix
    const auto xElements = 
        ranges::view::zip_with( std::multiplies< std::complex< double> >(),
                                row, xMultipliers );
    // the first row of the W matrix
    auto wElements = 
        xElements | ranges::view::transform(
                      [] ( const auto xValue ) -> std::complex< double >
                         { return std::complex< double >( 0., 2. ) * xValue; } );
    wElements.front() += 1.0;

    // the first row of the U matrix
    const auto uElements =
        ranges::view::zip_with( std::multiplies< std::complex< double> >(),
                                wElements, uMultipliers );

    // the exponential of the coulomb phase shift for the entrance channel
    const auto exponential = [&] {
      const auto wc = coulombShifts.front();
      return std::exp( std::complex< double >( 0.0, wc ) );
    }();

std::cout << "X(0,0) " << xElements[0] << std::endl;
std::cout << "W(0,0) " << wElements[0] << std::endl;
std::cout << "U(0,0) " << uElements[0] << std::endl;

    // the pi/k2 factor
    const auto factor = [&] {
      const auto waveNumber =
        this->incidentChannel().particlePair().waveNumber( energy );
      const auto squaredWaveNumber = waveNumber * waveNumber;
      const auto spinFactor = this->incidentChannel().statisticalSpinFactor();
      return pi / squaredWaveNumber * spinFactor;
    }();

// BEGIN REALLY BAD - GET TESTING GOING

    std::vector< Quantity< Barn > > result;
    Quantity< Barn > eliminated = 2. * factor * ( 1. - uElements.front().real() );
    for ( unsigned int i = 0; i < size; ++i ) {
      if ( i == 0 ) {
        result.push_back( factor * std::pow( std::abs( exponential - uElements[i] ), 2. ) );
      }
      else {
        result.push_back( factor * std::abs( std::pow( uElements[i], 2. ) ) );
      }
      eliminated -= result.back();
    }
    result.push_back( eliminated );
    return result;
// END REALLY BAD - GET TESTING GOING
  }
};
