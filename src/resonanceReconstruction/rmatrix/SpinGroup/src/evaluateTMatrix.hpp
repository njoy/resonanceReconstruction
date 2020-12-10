/**
 *  @brief Evaluate the elements of the T or X matrix at the given energy
 *
 *  The T or X matrix is defined as P^1/2 ( 1 - RL )^-1 R P^1/2 in which
 *  P is a diagonal matrix of the penetrabilities of each channel, R is the
 *  R matrix and L is a diagonal matrix defined as S - B + iP with S the shift
 *  factor and B the boundary condition of the channel.
 *
 *  @param[in] energy       the incident energy
 *  @param[in,out] result   a map containing the matrix elements
 */
void evaluateTMatrix(
         const Energy& energy,
         Map< ReactionChannelID, std::complex< double > >& result ) {

  // penetrability, sqrt(P) and channel identifiers for each channel
  const auto penetrabilities = this->penetrabilities( energy );
  const auto diagonalSqrtPMatrix = this->sqrtPenetrabilities( penetrabilities );
  const auto channels = this->channelIDs();

  // calculate the R_L = ( 1 - RL )^-1 R matrix
  decltype(auto) rlmatrix = this->rlmatrix_( energy,
                                             this->resonanceTable(),
                                             penetrabilities,
                                             this->channels() );

  // a lambda to process each channel
  const unsigned int size = channels.size();
  auto processChannel = [&] ( const unsigned int c ) {

    // the elements of the R_L = ( 1 - RL )^-1 R matrix for the current channel
    const auto row = ranges::make_iterator_range(
                        rlmatrix.data() + c * size,
                        rlmatrix.data() + ( c + 1 ) * size );

    // the row of the S or U matrix corresponding with the incident channel
    // S = U = Omega ( I + 2 i P^1/2 ( I - RL )^-1 R P^1/2 ) Omega
    // S = U = Omega ( I + 2 i P^1/2 R_L P^1/2 ) Omega
    // S = U = Omega ( I + 2 i T ) Omega
    // S = U = Omega W Omega
    const auto currentSqrtP = diagonalSqrtPMatrix[c];
    const auto elements =
        ranges::view::zip_with(
            [&] ( const auto tValue, const auto sqrtP )
                { return currentSqrtP * tValue * sqrtP; },
            row, diagonalSqrtPMatrix );

    // the element identifiers
    const auto current = channels[c];
    const auto identifiers =
        channels
          | ranges::view::transform(
                [&] ( const auto& id )
                    { return ReactionChannelID( current + "->" + id ); } );

    // assign into the map
    ranges::for_each(
      ranges::view::zip( identifiers, elements ),
      [&] ( const auto& pair ) -> void
          { result[ std::get< 0 >( pair ) ] = std::get< 1 >( pair ); } );
  };

  // process the channels
  const unsigned int start = 0;
  ranges::for_each( ranges::view::indices( start, size ), processChannel );
}
