ChannelRadii
makeChannelRadii( double ap, int naps, double awri, double neutronMass ) {

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
