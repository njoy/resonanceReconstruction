/**
 *  @brief Evaluate the elements of the T or X matrix at the given energy
 *
 *  @param[in] energy       the incident energy
 *  @param[in,out] result   a map containing the matrix elements
 */
void evaluateTMatrix(
         const Energy& energy,
         std::map< ReactionID, std::complex< double > >& result ) {

  ranges::for_each( this->groups_,
                    [&] ( auto& group )
                        { group.evaluateTMatrix( energy, result ); } );
}

