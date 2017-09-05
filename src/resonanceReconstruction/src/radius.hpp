inline auto radius( double scalar ){
  return
    [ r = scalar * rootBarn ]
    ( const Quantity<ElectronVolts> /* energy */ )
    { return r; };
}

inline auto radius( const ENDFtk::TAB1& tab1 ){
  using namespace interpolation;
  using namespace ranges;
  
  auto constructTable = []( auto&& region, int interpolant ){
    auto energy =
      region.first
      | view::transform( []( auto&& scalar ){ return scalar * electronVolts; } );

    auto radii =
      region.second
      | view::transform( []( auto&& scalar ){ return scalar * rootBarns; } );

    using Law1 =
      decltype( table::make< Histogram >( std::move(energy),
                                          std::move(radii) ) );

    using Law2 =
      decltype( table::make< LinearLinear >( std::move(energy),
                                             std::move(radii) ) );

    using Law3 =
      decltype( table::make< LinearLogarithmic >( std::move(energy),
                                                  std::move(radii) ) );

    using Law4 =
      decltype( table::make< LogarithmicLinear >( std::move(energy),
                                                  std::move(radii) ) );

    using Law5 =
      decltype( table::make< LogarithmicLogarithmic >( std::move(energy),
                                                       std::move(radii) ) );

    using ENDFvariant = Table< table::Variant< Law1, Law2, Law3, Law4, Law5 > >;

    switch( interpolant ){
    case 1: return ENDFvariant{ Law1{ std::move(energy), std::move(radii) } };
    case 2: return ENDFvariant{ Law2{ std::move(energy), std::move(radii) } };
    case 3: return ENDFvariant{ Law3{ std::move(energy), std::move(radii) } };
    case 4: return ENDFvariant{ Law4{ std::move(energy), std::move(radii) } };
    case 5: return ENDFvariant{ Law5{ std::move(energy), std::move(radii) } };
    }
    throw std::exception();
  };
  
  auto tables = 
    view::zip_with( constructTable, tab1.regions(), tab1.interpolants() )
    | to_vector;

  using ENDFvariant = decltype(tables)::value_type;
  
  auto table = Table< table::Vector< ENDFvariant >,
                      table::left::interval::Throws,
                      table::right::interval::Throws >{ std::move( tables ) };
  return table;
}
