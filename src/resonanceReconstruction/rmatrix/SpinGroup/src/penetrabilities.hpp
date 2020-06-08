auto penetrabilities( const Energy& energy ) const {
  
  return this->channels_
           | ranges::view::transform(
               [&] ( const auto& channel )
                   { return std::visit(
                       [&] ( const auto& channel )
                           { return channel.penetrability( energy ); },
                       channel ); } );
}
