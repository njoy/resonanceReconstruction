/**
 *  @brief Calculate possible values for a channel's total angular momentum J
 *
 *  The total angular momentum J for a channel can only have values  between
 *  abs(l - s) and l + s where l is the orbital momentum of the incoming wave
 *  and s is the channel spin (which in turn depends on the spin i of the
 *  incident particle and spin I of the target nucleus).
 *
 *  @param l   the orbital angular momentum of the incoming wave
 *  @param i   the spin of the incident particle
 *  @param I   the spin of the target nucleus
 *
 *  @return the possible values for the total angular momentum of the channel
 */
//template < typename Spin >
auto
possibleChannelTotalAngularMomentumValues( const OrbitalAngularMomentum& l,
                                           const Spin& i,
                                           const Spin& I ) {

  TotalAngularMomentum min = std::abs( std::abs( l - I ) - i );
  TotalAngularMomentum max = l + I + i;
  std::vector< TotalAngularMomentum > values = { min };
  while ( max > values.back() ) {

    values.push_back( values.back() + 1.0 );
  }
  return values;
}

/**
 *  @brief Calculate possible values for a channel's total angular momentum J
 *
 *  The total angular momentum J for a channel can only have values  between
 *  abs(l - s) and l + s where l is the orbital momentum of the incoming wave
 *  and s is the channel spin (which in turn depends on the spin i of the
 *  incident particle and spin I of the target nucleus).
 *
 *  @param l   the orbital angular momentum of the incoming wave
 *  @param s   the channel spin
 *
 *  @return the possible values for the total angular momentum of the channel
 */
auto
possibleChannelTotalAngularMomentumValues( const OrbitalAngularMomentum& l,
                                           const Spin& s ) {

  TotalAngularMomentum min = std::abs(l - s);
  TotalAngularMomentum max = l + s;
  std::vector< TotalAngularMomentum > values = { min };
  while ( max > values.back() ) {

    values.push_back( values.back() + 1.0 );
  }
  return values;
}
