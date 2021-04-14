auto phaseShifts( const Energy& energy ) const {

  const auto phaseShift = [=] ( const auto& channel ) {

    return channel.phaseShift( energy );
  };

  const auto getPhaseShift = [=] ( const auto& channel ) {

    return std::visit( phaseShift, channel );
  };

  return this->channels() | ranges::cpp20::views::transform( getPhaseShift );
}
