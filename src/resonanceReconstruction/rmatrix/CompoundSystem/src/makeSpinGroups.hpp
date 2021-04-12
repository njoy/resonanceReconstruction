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
  std::vector< std::pair< TotalAngularMomentum, Parity > > spins =
    ranges::to< std::vector< std::pair< TotalAngularMomentum, Parity > > >(
        channels | ranges::views::transform( getJpi ) )
    | ranges::actions::sort
    | ranges::actions::unique;

  // go over the Jpi values and create the spin groups
  for ( const auto& Jpi : spins ) {

    auto J = Jpi.first;
    auto parity = Jpi.second;

    auto filter = [&] ( const auto& channel ) {

      return ( channel.quantumNumbers().totalAngularMomentum() == J ) and
             ( channel.quantumNumbers().parity() == parity );
    };

    groups.emplace_back(
        ranges::to< std::vector< ParticleChannelData > >(
            channels | ranges::cpp20::views::filter( filter ) ) );
  }

  return groups;
}
