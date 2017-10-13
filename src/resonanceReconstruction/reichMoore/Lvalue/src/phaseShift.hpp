auto phaseShift( const double channelRatio ) const {
  return resonanceReconstruction::phaseShift
    ( this->orbitalAngularMomentum, channelRatio );
}
