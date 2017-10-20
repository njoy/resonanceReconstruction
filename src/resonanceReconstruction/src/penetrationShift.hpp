inline auto
penetrationShift( Integer<0> ){
  return []( const double channelRatio ) -> std::array< double, 2 > {
    const double p = ( channelRatio > 0 ) * channelRatio;
    const double s = 0.0;
    return {{ p, s }};
  };
}

inline auto
penetrationShift( Integer<1> ){
  return []( const double channelRatio ) -> std::array< double, 2 > {
    const bool positive = channelRatio > 0;
    const double rho2 = channelRatio * channelRatio;
    const double factor = 1. / ( 1 + rho2 );
    
    const double p = positive * factor * channelRatio * rho2;
    const double s = -factor;
    return {{ p, s }};
  };
}

inline auto
penetrationShift( Integer<2> ){
  return []( const double channelRatio ) -> std::array< double, 2 > {
    const bool positive = channelRatio > 0;
    const double rho2 = channelRatio * channelRatio;
    const double factor = 1. / ( 9. + rho2 * ( 3. + rho2 ) );
    
    const double p = positive * factor * channelRatio * rho2 * rho2;
    const double s = -factor * ( 18. + 3. * rho2 );
    return {{ p, s }};
  };
}

inline auto
penetrationShift( Integer<3> ){
  return []( const double channelRatio ) -> std::array< double, 2 > {
    const bool positive = channelRatio > 0;
    const double rho2 = channelRatio * channelRatio;
    const double factor = 1. / ( 225. + rho2 * ( 45. + rho2 * ( 6. + rho2 ) ) );
    
    const double p = positive * factor * channelRatio * rho2 * rho2 * rho2;
    const double s = -factor * ( 675. + rho2 * ( 90. + 6. * rho2 ) );
    return {{ p, s }};
  };
}

inline auto
penetrationShift( Integer<4> ){
  return []( const double channelRatio ) -> std::array< double, 2 > {
    const bool positive = channelRatio > 0;
    const double rho2 = channelRatio * channelRatio;
    const double rho4 = rho2 * rho2;
    const double factor = 1. / ( 11025. + rho2
                                 * ( 1575. + rho2
                                     * ( 135. + rho2
                                         * ( 10. + rho2 ) ) ) );

    const double p = positive * factor * channelRatio * rho4 * rho4;
    const double s = -factor * ( 44100. + rho2
                                 * ( 4725. + rho2
                                     * ( 270. + 10. * rho2 ) ) );
    return {{ p, s }};
  };
}

inline std::array< double, 2 >
penetrationShift( const int l, const double channelRatio ){
  const double rho2 = channelRatio * channelRatio;
  const double rho4 = rho2 * rho2;
  const bool positive = channelRatio > 0;
  
  // using horner's method
  switch(l){
  case 0: return {{ positive * channelRatio, 0.0 }};
  case 1: {
    const auto factor = 1. / ( 1 + rho2 );
    return {{ positive * factor * channelRatio * rho2, -factor }};
  }
  case 2: {
    const auto factor = 1. / ( 9. + rho2 * ( 3. + rho2 ) );
    return {{ positive * factor * channelRatio * rho4,
              -factor * ( 18. + 3. * rho2 ) }};
  }
  case 3: {
    const auto factor = 1. / ( 225. + rho2 * ( 45. + rho2 * ( 6. + rho2 ) ) );
    return {{ positive * factor * channelRatio * rho4 * rho2,
              -factor * ( 675. + rho2 * ( 90. + 6. * rho2 ) ) }};
  }
  case 4: {
    const auto factor = 1. / ( 11025. + rho2
                               * ( 1575. + rho2
                                   * ( 135. + rho2
                                       * ( 10. + rho2 ) ) ) );
    
    return {{ positive * factor * channelRatio * rho4 * rho4,
              -factor * ( 44100. + rho2
                          * ( 4725. + rho2
                              * ( 270. + 10. * rho2 ) ) ) }};
  }
  default: throw std::exception();
  }
}
