/**
 *  @brief Calculate possible values for the channel spin s
 *
 *  The channel spin s for a channel can only have values  between abs(i - I) 
 *  and i + I where i is the spin i of the incident particle (for a neutron that
 *  would be 0.5) and I is the spin of the target nucleus.
 *
 *  @param i   the spin of the incident particle
 *  @param I   the spin of the target nucleus
 *
 *  @return the possible values for the total angular momentum of the channel
 */
inline auto
possibleChannelSpinValues( const Spin& i, const Spin& I ) {

  Spin min = std::abs(i - I);
  Spin max = i + I;
  std::vector< Spin > values = { min };
  while ( max > values.back() ) {
    values.push_back( values.back() + 1.0 );
  }
  return values;
}

// use ranges::view::iota similar to ranges::view::closed_indices
