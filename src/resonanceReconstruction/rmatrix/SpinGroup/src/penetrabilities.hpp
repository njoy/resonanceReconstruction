auto penetrabilities( const Energy& energy ) const {

  auto penetrability = [=] ( const auto& channel )
                           { return channel.penetrability( energy ); };

  return this->channels()
           | ranges::view::transform(
                 [=] ( const auto& channel )
                     { return std::visit( penetrability, channel ); } );
}
