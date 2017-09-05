inline double phaseShift( const int l, const double channelRatio ){
  if ( not l ){ return channelRatio; }

  const auto rho2 = channelRatio * channelRatio;
  const auto offset = [&]{
    switch( l ){
    case 1: return channelRatio;
    case 2: return channelRatio * 3. / ( 3. - rho2 );
    case 3: return channelRatio * ( 15. - rho2 ) / ( 15. - 6. * rho2 );
    case 4: return channelRatio * ( 105. - 10. * rho2 )
                                  / ( 105. + rho2 * ( -45. + rho2 ) );
    default: throw std::exception();
    }
  }();
  
  return channelRatio - std::atan( offset );
}
