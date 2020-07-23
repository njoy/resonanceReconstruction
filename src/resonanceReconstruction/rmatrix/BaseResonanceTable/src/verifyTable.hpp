static
void verifyTable( const std::vector< ChannelID >& channels,
                  const std::vector< ResonanceType >& table ) {

  // verify that each resonance contains the same amount of resonances

  unsigned int expected = channels.size();

  const auto verifyNumberOfResonances = [expected] ( const auto& entry ) {

    unsigned int number = entry.widths().size();
    if ( expected != number ) {

      Log::error( "Number of reduced widths different from expected." );
      Log::info( "Found {}, expected {}", number, expected );
      throw std::exception();
    }
  };

  ranges::for_each( table, verifyNumberOfResonances );
}
