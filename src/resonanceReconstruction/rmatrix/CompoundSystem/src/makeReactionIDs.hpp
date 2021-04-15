static
std::vector< ReactionID >
makeReactionIDs( const std::vector< SpinGroup< Formalism, BoundaryOption > >& groups ) {

  std::vector< ReactionID > reactions;
  for ( const auto& group : groups ) {

    decltype(auto) groupreactions = group.reactionIDs();
    reactions.insert( reactions.end(),
                      groupreactions.begin(), groupreactions.end() );
  }

  ranges::cpp20::sort( reactions );
  reactions.erase( ranges::cpp20::unique( reactions ), reactions.end() );

  return reactions;
}
