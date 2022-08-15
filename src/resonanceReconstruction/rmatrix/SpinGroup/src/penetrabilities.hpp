auto penetrabilities( const Energy& energy ) const {

  const auto penetrability = [=] ( const auto& channel ) {

    return channel.penetrability( energy );
  };

  const auto getPenetrability = [=] ( const auto& channel ) {

    return std::visit( penetrability, channel );
  };

  return this->channels() | ranges::cpp20::views::transform( getPenetrability );
}
