std::optional< ChannelRadiusTable >
makeChannelRadiusTable( const std::optional< ENDF::ScatteringRadius >& radius ) {

  if ( radius ) {

    auto makeTable = [] ( auto&& region, int interpolant )
      -> TableVariant< Energy, ChannelRadius > {

      auto toEnergy = [] ( const auto& value ) { return value * electronVolt; };
      auto toRadius = [] ( const auto& value ) { return value * rootBarn; };

      std::vector< Energy > energies =
          region.first | ranges::view::transform( toEnergy );
      std::vector< ChannelRadius > radii =
          region.second | ranges::view::transform( toRadius );

      switch( interpolant ) {

        case 1: {

          return HistogramTable< Energy, ChannelRadius >( std::move( energies ),
                                                          std::move( radii ) );
        }
        case 2: {

          return LinLinTable< Energy, ChannelRadius >( std::move( energies ),
                                                       std::move( radii ) );
        }
        case 3: {

          return LinLogTable< Energy, ChannelRadius >( std::move( energies ),
                                                       std::move( radii ) );
        }
        case 4: {

          return LogLinTable< Energy, ChannelRadius >( std::move( energies ),
                                                       std::move( radii ) );
        }
        case 5: {

          return LogLogTable< Energy, ChannelRadius >( std::move( energies ),
                                                       std::move( radii ) );
        }
        default : {

          throw std::runtime_error( "You somehow reached unreachable code" );
        }
      }
    };

    const auto regions = radius->regions();
    const auto interpolants = radius->interpolants();
    std::vector< TableVariant< Energy, ChannelRadius > > tables =
      ranges::view::zip_with( makeTable, regions, interpolants );

    ChannelRadiusTable table(
        MultiRegionTable< Energy, ChannelRadius >( std::move( tables ) ) );
    return std::make_optional( std::move( table ) );
  }
  return std::nullopt;
}
