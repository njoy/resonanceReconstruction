auto coulombShifts( const Energy& energy ) const {

  const auto coulombPhaseShift = [=] ( const auto& channel ) {

    return channel.coulombPhaseShift( energy );
  };

  return this->channels()
           | ranges::cpp20::views::transform(
                 [=] ( const auto& channel )
                     { return std::visit( coulombPhaseShift, channel ); } );
}
