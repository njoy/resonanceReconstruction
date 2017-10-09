auto penetrationShift( const double channelRatio ) const {
  return resonanceReconstruction::penetrationShift
    ( this->orbitalAngularMomentum, channelRatio );
}
