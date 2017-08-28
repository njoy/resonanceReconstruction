auto
penetrationShift( const double channelRatio ) const {
  return singleLevelBreitWigner::penetrationShift( this->orbitalAngularMomentum,
                                                   channelRatio );
}
