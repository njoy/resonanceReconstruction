static
std::vector< SpinGroup< Formalism, BoundaryOption > >
makeSpinGroups( std::vector< ParticleChannelData >&& channels ) {

  std::vector< SpinGroup< Formalism, BoundaryOption > > groups;

  // usefull lambdas
  auto getJpi = [] ( const auto& channel ) {

    return std::make_pair( channel.quantumNumbers().totalAngularMomentum(),
                           channel.quantumNumbers().parity() );
  };

  // get the different Jpi values in the channels
  using Pair = std::pair< TotalAngularMomentum, Parity >;
  std::vector< Pair > spins =
      ranges::to< std::vector< Pair > >(
          channels | ranges::cpp20::views::transform( getJpi ) );

  ranges::cpp20::sort( spins );
  spins.erase( ranges::cpp20::unique( spins ), spins.end() );

  // go over the Jpi values and create the spin groups
  for ( const auto& Jpi : spins ) {

    auto J = Jpi.first;
    auto parity = Jpi.second;

    auto filter = [&] ( const auto& channel ) {

      return ( channel.quantumNumbers().totalAngularMomentum() == J ) and
             ( channel.quantumNumbers().parity() == parity );
    };

    auto filtered = channels | ranges::cpp20::views::filter( filter );
    groups.emplace_back(
        std::vector< ParticleChannelData >( ranges::cpp20::begin( filtered ),
                                            ranges::cpp20::end( filtered ) ) );
  }

  return groups;
}
