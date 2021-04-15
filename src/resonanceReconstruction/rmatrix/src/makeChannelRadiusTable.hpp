std::optional< ChannelRadiusTable >
makeChannelRadiusTable( const std::optional< endf::ScatteringRadius >& radius ) {

  if ( radius ) {

    auto makeTable = [] ( auto&& region, int interpolant )
      -> TableVariant< Energy, ChannelRadius > {

      const auto toEnergy = [] ( const auto& value ) -> Energy
                               { return value * electronVolt; };
      const auto toRadius = [] ( const auto& value ) -> ChannelRadius
                               { return value * rootBarn; };

      std::vector< Energy > energies =
          ranges::to< std::vector< Energy > >(
              region.first | ranges::cpp20::views::transform( toEnergy ) );
      std::vector< ChannelRadius > radii =
          ranges::to< std::vector< ChannelRadius > >(
              region.second | ranges::cpp20::views::transform( toRadius ) );

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
      ranges::to< std::vector< TableVariant< Energy, ChannelRadius > > >(
        ranges::views::zip_with( makeTable, regions, interpolants ) );

    ChannelRadiusTable table(
        MultiRegionTable< Energy, ChannelRadius >( std::move( tables ) ) );
    return std::make_optional( std::move( table ) );
  }
  return std::nullopt;
}
