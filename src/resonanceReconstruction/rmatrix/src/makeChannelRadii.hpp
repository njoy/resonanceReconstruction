ChannelRadii
makeChannelRadii( double ap, const std::optional< ChannelRadiusTable >& nro,
                  int naps, double awri, double neutronMass ) {

  if ( nro ) {

    switch ( naps ) {

      case 0 : {

        return ChannelRadii(
                 ( 0.123 * std::pow( awri * neutronMass, 1. / 3. )
                   + 0.08 ) * rootBarn,
                 nro.value() );
      }
      case 1 : {

        return ChannelRadii( nro.value() );
      }
      case 2 : {

        return ChannelRadii( ap * rootBarn, nro.value() );
      }
      default : {

        throw std::runtime_error( "You somehow reached unreachable code" );
      }
    }
  }
  else {

    if ( not naps ) {

      return ChannelRadii(
               ( 0.123 * std::pow( awri * neutronMass, 1. / 3. )
                 + 0.08 ) * rootBarn,
               ap * rootBarn );
    }
    else {

      return ChannelRadii( ap * rootBarn );
    }
  }
}
