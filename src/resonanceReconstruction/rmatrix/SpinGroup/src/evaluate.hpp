/**
 *  @brief Evaluate the cross sections at the given energy
 *
 *  @param[in] energy       the incident energy
 *  @param[in,out] result   a map containing the accumulated cross sections
 */
void evaluate( const Energy& energy,
               Map< ReactionID, CrossSection >& result ) {

  // penetrability, Coulomb phase shift, sqrt(P) and Omega = exp( i(w - phi) )
  // for each channel except the eliminated capture channel
  const auto penetrabilities = this->penetrabilities( energy );
  const auto coulombShifts = this->coulombShifts( energy );
  const auto diagonalSqrtPMatrix = this->sqrtPenetrabilities( penetrabilities );
  const auto diagonalOmegaMatrix = this->omegas( energy, coulombShifts );

  // calculate the R_L = ( 1 - RL )^-1 R matrix
  decltype(auto) rlmatrix = this->rlmatrix_( energy,
                                             this->resonanceTable(),
                                             penetrabilities,
                                             this->channels() );

  // the pi/k2 * gJ factor
  const auto factor = [&] {
    auto factor = [&] ( const auto& channel ) {
      const auto waveNumber = channel.waveNumber( energy );
      const auto squaredWaveNumber = waveNumber * waveNumber;
      const auto spinFactor = channel.statisticalSpinFactor();
      return pi / squaredWaveNumber * spinFactor;
    };
    return std::visit( factor, this->channels_[ this->incident_.front() ] );
  }();

  // the cross section identifiers
  const auto identifiers = this->reactionIDs();

  // a lambda to process each incident channel
  auto processIncidentChannel = [&] ( const unsigned int c ) {

    // lambda to derive a kronecker delta array for the current incident channel
    const unsigned int size = this->channels().size();
    auto delta = [c,size] ( const auto value ) {
      return ranges::views::concat(
                 ranges::views::repeat_n( 0., c ),
                 ranges::cpp20::views::single( value ),
                 ranges::views::repeat_n( 0., size - c - 1 ) );
    };

    // the elements of the R_L = ( 1 - RL )^-1 R matrix for the incident channel
    const auto row = ranges::make_subrange(
                        rlmatrix.data() + c * size,
                        rlmatrix.data() + ( c + 1 ) * size );

    // the row of the S or U matrix corresponding with the incident channel
    // S = U = Omega ( I + 2 i P^1/2 ( I - RL )^-1 R P^1/2 ) Omega
    // S = U = Omega ( I + 2 i P^1/2 R_L P^1/2 ) Omega
    // S = U = Omega ( I + 2 i T ) Omega
    // S = U = Omega W Omega
    const auto incidentSqrtP = diagonalSqrtPMatrix[c];
    const auto incidentOmega = diagonalOmegaMatrix[c];
    const auto uElements =
        ranges::views::zip_with(
            [&] ( const auto delta, const auto tValue,
                  const auto sqrtP, const auto omega )
                { return incidentOmega *
                         ( delta + std::complex< double >( 0., 2. ) *
                                   incidentSqrtP * tValue * sqrtP ) * omega; },
            delta( 1.0 ), row, diagonalSqrtPMatrix, diagonalOmegaMatrix );

    // the exponential of the coulomb phase shift for the incident channel
    const auto exponential =
      std::exp( std::complex< double >( 0., coulombShifts[c] ) );

    // the cross section values for channel c to c' - independent of formalism
    // sigma_cc' = norm( exp( iw_c ) delta_cc' - U_cc' )
    const auto sigma =
      ranges::views::zip_with(
          [&] ( const auto delta, const auto uValue )
              { return std::norm( delta - uValue ); },
          delta( exponential ),
          uElements );

    // the eliminated capture channel - Reich-Moore only
    const auto capture =
      ranges::cpp20::views::single(
          ranges::accumulate(
              uElements | ranges::cpp20::views::transform(
                              [] ( const auto value ) -> double
                                 { return std::norm( value ); } ),
              1., ranges::minus() ) );

    // concat and multiply by pi / k^2 g_J
    const auto crossSections =
      ranges::views::concat( sigma, capture )
        | ranges::cpp20::views::transform(
              [=] ( const auto value ) -> Quantity< Barn >
                  { return factor * value; } );

    // accumulate results in the map
    ranges::cpp20::for_each(
      ranges::views::zip( identifiers, crossSections ),
      [&] ( const auto& pair ) -> void
          { result[ std::get< 0 >( pair ) ] += std::get< 1 >( pair ); } );
  };

  // process the incident channels
  ranges::cpp20::for_each( this->incident_, processIncidentChannel );
}
