auto phaseShifts( const Energy& energy ) const {

  auto phaseShift = [=] ( const auto& channel )
                        { return channel.phaseShift( energy ); };

  return this->channels()
           | ranges::views::transform(
                 [=] ( const auto& channel )
                     { return std::visit( phaseShift, channel ); } );
}
