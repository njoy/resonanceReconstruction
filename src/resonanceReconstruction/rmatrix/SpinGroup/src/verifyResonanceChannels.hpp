static
void verifyResonanceChannels( const std::vector< ParticleChannel >& channels,
                              const ResonanceTable& table ) {

  // check size
  const auto size = channels.size();
  const auto number = table.numberChannels();
  if ( size != number ) {

    Log::error( "Inconsistent size of resonance table with respect to the "
                "channels in the spin group." );
    Log::info( "Expected {}, found {}", size, number );
    throw std::exception();
  }

  // verify that the channel labels are in the same order

  const auto getChannelID = [] ( const auto& entry ) {

    return std::visit( [] ( const auto& channel )
                          { return channel.channelID(); },
                       entry );
  };

  const auto checkLabels = [&] ( const auto& pair ) {

    const auto channel = std::get< 0 >( pair );
    const auto resonance = std::get< 1 >( pair );
    if ( channel != resonance ) {

      Log::error( "Channels and resonance table are not in the same order." );
      Log::info( "Expected channel {}, found {} in the resonance table",
                 channel, resonance );
      throw std::exception();
    }
  };

  ranges::cpp20::for_each(
      ranges::views::zip( channels | ranges::views::transform( getChannelID ),
                          table.channels() ),
      checkLabels );
}
