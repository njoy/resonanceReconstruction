std::array< double, 2 >
penetrationShift( const double channelRatio ) const {
  const double rho2 = channelRatio * channelRatio;
  const double rho4 = rho2 * rho2;

  // using horner's method
  auto pair = [&]() -> std::array< double, 2 > {
    switch( this->orbitalAngularMomentum ){
    case 0: return {{ channelRatio, 0.0 }};
    case 1: {
      const auto factor = 1. / ( 1 + rho2 ) ;
      return {{ factor * channelRatio * rho2, -factor }};
    } case 2: {
      const auto factor = 1. / ( 9. + rho2 * ( 3. + rho2 ) );
      return {{ factor * channelRatio * rho4, -factor * ( 18. + 3. * rho2 ) }};
    } case 3: {
      const auto factor = 1. / ( 225. + rho2 * ( 45. + rho2 * ( 6. + rho2 ) ) );
      return {{ factor * channelRatio * rho4 * rho2,
                -factor * ( 675. + rho2 * ( 90. + 6. * rho2 ) ) }};
    } case 4: {
      const auto factor =
        1. / ( 11025. + rho2
                        * ( 1575. + rho2
                                    * ( 135. + rho2
                                               * ( 10. + rho2 ) ) ) );
      return
      {{ factor * channelRatio * rho4 * rho4,
         -factor * ( 44100. + rho2 * ( 4725. + rho2 * ( 270. + 10. * rho2 ) ) ) }};
    } default: throw std::exception();
    }
  }();
  
  if ( channelRatio < 0 ){ pair[0] = 0.0; }
  return pair;
}
