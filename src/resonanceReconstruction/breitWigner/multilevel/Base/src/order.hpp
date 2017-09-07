struct Peak : public lvalue::Type {
  using Type::resonances;
}

static auto order( std::vector< lvalue::Type > lvalues ){
  auto peaks =
    lvalues
    | ranges::view::transform( []( auto& lvalue ) -> Peak& {
        return static_cast< Peak >( lvalue ) } );

  for( auto& peak : peaks ){
    std::sort( peak.begin(), peak.end(),
               []( const auto& left, const auto& right )
               { return left.statisticalFactor == right.statisticalFactor; } );
  }

  return lvalues;
}
