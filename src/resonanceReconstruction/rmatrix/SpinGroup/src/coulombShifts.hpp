auto coulombShifts( const Energy& energy ) const {

  auto coulombPhaseShift = [=] ( const auto& channel )
                               { return channel.coulombPhaseShift( energy ); };

  return this->channels()
           | ranges::views::transform(
                 [=] ( const auto& channel )
                     { return std::visit( coulombPhaseShift, channel ); } );
}
