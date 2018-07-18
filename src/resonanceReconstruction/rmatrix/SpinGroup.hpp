/**
 *  @class
 *  @brief A spin group corresponding to a J,pi value
 */
class SpinGroup {

  /* fields */
  std::vector< Channel > channels_; // first channel is entrance
  ResonanceTable parameters_;

// BEGIN REALLY BAD - GET TESTING GOING
  void matrices( const Energy& energy,
                 Matrix< double >& pSqrtMatrix,
                 Matrix< std::complex< double > >& lMatrix,
                 Matrix< std::complex< double > >& omegaMatrix ) const {

    const unsigned int size = this->parameters_.numberChannels();
    pSqrtMatrix = Matrix< double >::Zero( size, size );
    lMatrix = Matrix< std::complex< double > >::Zero( size, size );
    omegaMatrix = Matrix< std::complex< double > >::Zero( size, size );
    for ( unsigned int i = 0; i < size; ++i ) {

      const auto P = this->channels_[i].penetrability( energy );
      const auto S = this->channels_[i].shiftFactor( energy );
      const auto B = this->channels_[i].boundaryCondition();
      const auto w = this->channels_[i].coulombPhaseShift( energy );
      const auto phi = this->channels_[i].phaseShift( energy );

      pSqrtMatrix.diagonal()[i] = std::sqrt( P );
      lMatrix.diagonal()[i] = std::complex< double >( S - B, P );
      omegaMatrix.diagonal()[i] = std::exp( std::complex< double >( 0.0, w - phi ) );
    }
  }
// END REALLY BAD - GET TESTING GOING

public:

  /* constructor */
  SpinGroup( std::vector< Channel >&& channels, ResonanceTable&& table ) :
    channels_( std::move( channels ) ), parameters_( std::move( table ) ) {}

  auto resonanceTable() const { return this->parameters_; }
  auto channels() const { return ranges::view::all( this->channels_ ); }
  auto incidentChannel() const { return this->channels_.front(); }
  auto outgoingChannels() const {
    return ranges::view::drop_exactly( this->channels(), 1 );
  }

  auto evaluate( const Energy& energy ) const {

    const unsigned int size = this->parameters_.numberChannels();

    // get the R-matrix
    Matrix< std::complex< double > > rMatrix =
        this->resonanceTable().rmatrix( energy );

// BEGIN REALLY BAD - GET TESTING GOING

    //! @todo we only need the first row of ( L^-1 - R )^-1 and each column of R

    // get the other matrices
    Matrix< double > pSqrtMatrix( size, size );
    Matrix< std::complex< double > > lMatrix( size, size );
    Matrix< std::complex< double > > omegaMatrix( size, size );
    this->matrices( energy, pSqrtMatrix, lMatrix, omegaMatrix );

    // calculate the U matrix
    Matrix< std::complex< double > > uMatrix =
        omegaMatrix *
          ( Matrix< double >::Identity( size, size ) + 
            std::complex< double >( 0.0, 2.0 ) * pSqrtMatrix *
            ( Matrix< double >::Identity( size, size )
              - rMatrix * lMatrix ).inverse() * rMatrix * pSqrtMatrix ) *
        omegaMatrix;

// END REALLY BAD - GET TESTING GOING

    // get the pi/k2 factor
    const auto factor = [&] {
      const auto waveNumber =
        this->incidentChannel().particlePair().waveNumber( energy );
      const auto squaredWaveNumber = waveNumber * waveNumber;
      const auto spinFactor = this->incidentChannel().statisticalSpinFactor();
      return pi / squaredWaveNumber * spinFactor;
    }();

// BEGIN REALLY BAD - GET TESTING GOING
    std::vector< Quantity< Barn > > result;
    Quantity< Barn > eliminated = 2. * factor * ( 1. - uMatrix( 0, 0 ).real() );
    for ( unsigned int i = 0; i < size; ++i ) {
      if ( i == 0 ) {
        const double wc = this->channels_[i].coulombPhaseShift( energy );
        result.push_back( factor * std::pow( std::abs( std::exp( std::complex< double >( 0.0, wc ) ) - uMatrix( i, i ) ), 2. ) );
      }
      else {
        result.push_back( factor * std::abs( std::pow( uMatrix( i, i ), 2. ) ) );
      }
      eliminated -= result.back();
    }
    result.push_back( eliminated );
    return result;
// END REALLY BAD - GET TESTING GOING
  }
};
