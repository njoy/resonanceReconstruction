#ifndef NJOY_R2_PYTHON_CONVERSION
#define NJOY_R2_PYTHON_CONVERSION

// system includes

// other includes
#include "resonanceReconstruction/Quantity.hpp"

template < typename Unit, typename Magnitude >
using Quantity = dimwits::Quantity< Unit, Magnitude >;
using AtomicMass = njoy::resonanceReconstruction::AtomicMass;
using ElectricalCharge = njoy::resonanceReconstruction::ElectricalCharge;

template < typename Unit, typename Magnitude >
Magnitude removeUnit( const Quantity< Unit, Magnitude >& quantity ) {

  return quantity.value;
}

inline AtomicMass toAtomicMass( double value ) {

  return value * dimwits::daltons;
}

inline ElectricalCharge toElectricalCharge( double value ) {

  return value * dimwits::coulombs;
}

#endif
