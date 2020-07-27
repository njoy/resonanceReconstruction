/**
 *  @brief Evaluate the cross sections at the given energy
 *
 *  @param[in] energy       the incident energy
 *  @param[in,out] result   a map containing the accumulated cross sections
 */
void evaluate( const Energy& energy,
               std::map< ReactionID, CrossSection >& result ) {

  ranges::for_each( this->groups_,
                    [&] ( auto& group )
                        { group.evaluate( energy, result ); } );
}
