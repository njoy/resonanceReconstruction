static
std::vector< ReactionID >
makeReactionIDs( const std::vector< SpinGroupType >& groups ) {

  std::vector< ReactionID > reactions;

  reactions = groups | ranges::view::transform(
                           [] ( const auto& group )
                              { return group.reactionIDs(); } )
                     | ranges::view::join;

  std::sort( reactions.begin(), reactions.end() );
  reactions.erase( std::unique( reactions.begin(), reactions.end() ),
                   reactions.end() );

  return reactions;
}
