inline double phaseShift( const int L, const double channelRatio ){
  if ( L == 0 )
    return channelRatio;

  const double rho2 = channelRatio * channelRatio;
  double offset;
  switch( L ){
    case 1: offset = 1; break;
    case 2: offset = 3/(3-rho2); break;
    case 3: offset = (15-rho2)/(15-6*rho2); break;
    case 4: offset = (105-10*rho2)/(105+rho2*(rho2-45)); break;
    default: throw std::exception();
  };
  
  return channelRatio - std::atan( channelRatio * offset );
}
