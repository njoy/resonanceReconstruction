/**
 *  @brief Evaluate the elements of the T or X matrix at the given energy
 *
 *  @param[in] energy       the incident energy
 *  @param[in,out] result   a map containing the matrix elements
 */
void evaluateTMatrix(
         const Energy& energy,
         tsl::hopscotch_map< ReactionID, Quantity< Barn > >& result ) {

  // penetrability, sqrt(P) and channel identifiers for each channel
  const auto penetrabilities = this->penetrabilities( energy );
  const auto diagonalSqrtPMatrix = this->sqrtPenetrabilities( penetrabilities );
  const auto channels = this->channelIDs();

  // calculate the diagonal of the L matrix based on the BoundaryOption template
  this->diagonalLMatrix_ =
      this->calculateLDiagonal( energy, penetrabilities, BoundaryOption() );

  // calculate the R_L = ( 1 - RL )^-1 R matrix based on the Formalism template
  this->calculateRLMatrix( energy, Formalism() );

  // a lambda to process each channel
  const unsigned int size = channels.size();
  auto processChannel = [&] ( const unsigned int c ) {

    // the elements of the ( 1 - RL )^-1 R matrix for the current channel
    const auto row =
    ranges::make_iterator_range(
        this->matrix_.row(c).data(),
        this->matrix_.row(c).data() + size );

    // the row of the T or X matrix corresponding with the current channel
    const auto currentSqrtP = diagonalSqrtPMatrix[c];
    const auto elements =
        ranges::view::zip_with( 
            [&] ( const auto tValue, const auto sqrtP )
                { return currentSqrtP * tValue * sqrtP; },
            row, diagonalSqrtPMatrix );

    // the element identifiers
    const auto current = diagonalSqrtPMatrix[c];
    const auto identifiers =
        channels
          | ranges::view::transform(
                [&] ( const auto& id )
                    { return ReactionID( current + "->" + id ); } );

    // accumulate results
    ranges::for_each(
      ranges::view::zip( identifiers, elements ),
      [&] ( const auto& pair ) -> void
          { result[ std::get< 0 >( pair ) ] = std::get< 1 >( pair ); } );
  };

  // process the channels
  const unsigned int start = 0;
  ranges::for_each( ranges::view::indices( start, size ), processChannel );
}

