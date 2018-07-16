/**
 *  @class
 *  @brief A spin group corresponding to a J,pi value
 */
class SpinGroup {

  /* fields */
  std::vector< Channel > channels_;
  ResonanceTable parameters_;

  void matrices( const Energy& energy,
                 Matrix< std::complex< double > > pSqrtMatrix,
                 Matrix< std::complex< double > > lInvMatrix,
                 Matrix< std::complex< double > > omegaMatrix ) const {

    const unsigned int size = this->parameters_.numberChannels();
    pSqrtMatrix = Matrix< std::complex< double > >::Zero( size, size );
    lInvMatrix = Matrix< std::complex< double > >::Zero( size, size );
    omegaMatrix = Matrix< std::complex< double > >::Zero( size, size );
    for ( unsigned int i = 0; i < size; ++i ) {

      const auto P = this->channels_[i].penetrability( energy );
      const auto S = this->channels_[i].shiftFactor( energy );
      const auto B = this->channels_[i].boundaryCondition();
      const auto w = this->channels_[i].coulombPhaseShift( energy );
      const auto phi = this->channels_[i].phaseShift( energy );

      pSqrtMatrix.diagonal()[i] = std::sqrt( P );
      lInvMatrix.diagonal()[i] = 1. / std::complex< double >( S - B, P );
      omegaMatrix.diagonal()[i] = std::exp( std::complex< double >( 0.0, w - phi ) );
    }
  }

public:

  /* constructor */

  auto resonanceTable() const { return this->parameters_; }
  auto channels() const { return ranges::view::all( this->channels_ ); }

  Matrix< std::complex< double > > generateCollisionMatrix( const Energy& energy ) const {

    const unsigned int size = this->parameters_.numberChannels();

    // get the R-matrix
    Matrix< std::complex< double > > rMatrix =
        this->resonanceTable().rmatrix( energy );

    // get the other matrices
    Matrix< std::complex< double > > pSqrtMatrix( size, size );
    Matrix< std::complex< double > > lInvMatrix( size, size );
    Matrix< std::complex< double > > omegaMatrix( size, size );
    this->matrices( energy, pSqrtMatrix, lInvMatrix, omegaMatrix );

    // calculate the U matrix
    return
        omegaMatrix * ( Matrix< double >::Identity( size, size ) + 
                        std::complex< double >( 0.0, 2.0 ) * pSqrtMatrix *
                        lInvMatrix * ( lInvMatrix - rMatrix ).inverse() *
                        rMatrix * pSqrtMatrix ) * omegaMatrix;
  }
};
