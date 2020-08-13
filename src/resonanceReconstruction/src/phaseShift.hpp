inline double phaseShift( const int L, const double channelRatio ){
  if ( L == 0 )
    return channelRatio;

  const double rho2 = channelRatio * channelRatio;
  const double offset = [&]{
    switch( L ){
      case 1: return 1.0;
      case 2: return 3/(3-rho2);
      case 3: return (15-rho2)/(15-6*rho2);
      case 4: return (105-10*rho2)/(105+rho2*(rho2-45));
      default: throw std::exception();
    }
  }();
  
  return channelRatio - std::atan( channelRatio * offset );
}
