static
ChannelID makeChannelID( const std::string& id,
                         const ChannelQuantumNumbers& numbers ) {

  return id + numbers.toString();
}
