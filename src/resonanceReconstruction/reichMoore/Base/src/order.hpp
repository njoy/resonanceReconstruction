struct Peak : public Lvalue {
  using Lvalue::resonances;
};

static auto order( std::vector< Lvalue > lvalues ){
  auto peaks =
    lvalues
    | ranges::views::transform
      ( []( auto& lvalue ) -> std::vector< Resonance >&
        { return static_cast< Peak& >( lvalue ).resonances; } );

  for( auto& peak : peaks ){
    std::sort( peak.begin(), peak.end(),
               []( const auto& left, const auto& right )
               { return left.flaggedStatisticalFactor
                        < right.flaggedStatisticalFactor; } );
  }

  return lvalues;
}
