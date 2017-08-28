inline double phaseShift( const int l, const double channelRatio ){
  const auto rho2 = channelRatio * channelRatio;
  const auto offset = [&]{
    switch( l ){
    case 0: return 0.0;
    case 1: return std::atan( channelRatio );
    case 2: return std::atan( 3. * channelRatio / ( 3. - rho2 ) );
    case 3: {
      return std::atan( channelRatio * ( 15. - rho2 ) / ( 15. - 6. * rho2 ) );
    }
    case 4: {
      return std::atan( channelRatio * ( 105. - 10. * rho2 )
                        / ( 105. + rho2 * ( -45. + rho2 ) ) );

    }
    default: throw std::exception();
    }
  }();
  return channelRatio - offset;
}
