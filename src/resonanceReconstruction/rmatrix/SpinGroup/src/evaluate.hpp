/**
 *  @brief Evaluate the cross sections at the given energy
 *
 *  @param[in] energy       the incident energy
 *  @param[in,out] result   a map containing the accumulated cross sections
 */
void evaluate( const Energy& energy,
               tsl::hopscotch_map< ReactionID, Quantity< Barn > >& result ) {

  // penetrability, shift factor, phase shift and Coulomb phase shift
  // for each channel except the eliminated capture channel
  const auto penetrabilities = this->penetrabilities( energy );
  const auto coulombShifts = this->coulombShifts( energy );
  const auto phaseShifts = this->phaseShifts( energy );

  // calculate the diagonal of the L matrix based on the BoundaryOption template
  this->diagonalLMatrix_ =
      this->calculateLDiagonal( energy, penetrabilities, BoundaryOption() );

  // the diagonal of the sqrt(P) matrix
  const auto diagonalSqrtPMatrix =
    penetrabilities
      | ranges::view::transform(
            [] ( const double penetrability ) -> double
               { return std::sqrt( penetrability ); } );

  // the diagonal of the Omega matrix = exp( i ( w - phi ) )
  const auto diagonalOmegaMatrix =
    ranges::view::zip_with(
        [] ( const double w, const double phi )
           { return std::exp( std::complex< double >( 0.0, w - phi ) ); },
        coulombShifts, phaseShifts );

  // calculate the R_L = ( 1 - RL )^-1 R matrix based on the Formalism template
  this->calculateRLMatrix( energy, Formalism() );

  // the pi/k2 * gJ factor
  const auto factor = [&] {
    auto factor = [&] ( const auto& channel ) {
      const auto waveNumber = channel.particlePair().waveNumber( energy );
      const auto squaredWaveNumber = waveNumber * waveNumber;
      const auto spinFactor = channel.statisticalSpinFactor();
      return pi / squaredWaveNumber * spinFactor;
    };
    return std::visit( factor, this->channels_[ this->incident_.front() ] );
  }();

  // the cross section identifiers
  const auto identifiers = this->reactions();

  // a lambda to process each incident channel
  auto processIncidentChannel = [&] ( const unsigned int c ) {

    // lambda to derive a kronecker delta array for the current incident channel
    const unsigned int size = this->channels().size();
    auto delta = [c,size] ( const auto value ) {
      return ranges::view::concat(
                 ranges::view::repeat_n( 0., c ),
                 ranges::view::single( value ),
                 ranges::view::repeat_n( 0., size - c - 1 ) );
    };

    // the elements of the ( 1 - RL )^-1 R matrix for each channel (assumes
    // the incident channel is the first channel)
    const auto row =
    ranges::make_iterator_range(
        this->matrix_.row(c).data(),
        this->matrix_.row(c).data() + size );

    // the row of the U matrix corresponding with the incident channel
    const auto incidentSqrtP = diagonalSqrtPMatrix[c];
    const auto incidentOmega = diagonalOmegaMatrix[c];
    const auto uElements =
      ranges::view::zip_with( 
          [&] ( const auto delta, const auto tValue,
                const auto sqrtP, const auto omega )
              { return incidentOmega *
                       ( delta + std::complex< double >( 0., 2. ) *
                                 incidentSqrtP * tValue * sqrtP ) * omega; },
          delta( 1.0 ), row, diagonalSqrtPMatrix, diagonalOmegaMatrix );

    // the exponential of the coulomb phase shift for the incident channel
    const auto exponential =
      std::exp( std::complex< double >( 0., coulombShifts[c] ) );

    // lambda to calculate a norm squared
    auto normSquared = [] ( const auto value ) -> double
                          { return std::pow( std::abs( value ), 2. ); };

    // the cross section values
    const auto crossSections =
      ranges::view::concat(
          ranges::view::zip_with(
              [&] ( const auto delta, const auto uValue )
                  { return normSquared( delta - uValue ); },
              delta( exponential ),
              uElements ),
          ranges::view::single(
              ranges::accumulate(
                  uElements | ranges::view::transform( normSquared ), 1.,
                  ranges::minus() ) ) )
        | ranges::view::transform(
              [=] ( const auto value ) -> Quantity< Barn >
                  { return factor * value; } );

    // accumulate results
    //! @todo debug build requires accumulating cross section values inside the
    //        processIncidentChannel lambda
    ranges::for_each(
      ranges::view::zip( identifiers, crossSections ),
      [&] ( const auto& pair ) -> void
          { result[ std::get< 0 >( pair ) ] += std::get< 1 >( pair ); } );
  };

  // process the incident channels
  ranges::for_each( this->incident_, processIncidentChannel );
}

