auto coulombShifts( const Energy& energy ) const {

  auto coulombPhaseShift = [=] ( const auto& channel )
                               { return channel.coulombPhaseShift( energy ); };

  return this->channels()
           | ranges::view::transform(
                 [=] ( const auto& channel )
                     { return std::visit( coulombPhaseShift, channel ); } );
}
