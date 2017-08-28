double phaseShift( const double channelRatio ) const {
  return singleLevelBreitWigner::phaseShift( this->orbitalAngularMomentum,
                                             channelRatio );
}
