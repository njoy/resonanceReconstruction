#ifndef NJOY_R2_PYTHON_CONVERSION
#define NJOY_R2_PYTHON_CONVERSION

// system includes

// other includes
#include "resonanceReconstruction/Quantity.hpp"

template < typename Unit, typename Magnitude >
using Quantity = dimwits::Quantity< Unit, Magnitude >;
using AtomicMass = njoy::resonanceReconstruction::AtomicMass;
using ElectricalCharge = njoy::resonanceReconstruction::ElectricalCharge;
using ChannelRadius = njoy::resonanceReconstruction::ChannelRadius;
using Energy = njoy::resonanceReconstruction::Energy;
using ReducedWidth = njoy::resonanceReconstruction::ReducedWidth;
using QValue = njoy::resonanceReconstruction::QValue;

template < typename Unit, typename Magnitude >
Magnitude removeUnit( const Quantity< Unit, Magnitude >& quantity ) {

  return quantity.value;
}

template < typename Unit, typename Magnitude >
std::vector< double > removeArrayUnit( const std::vector< Quantity< Unit, Magnitude > >& array ) {

  return quantity.value;
}

inline AtomicMass toAtomicMass( double value ) {

  return value * dimwits::daltons;
}

inline ElectricalCharge toElectricalCharge( double value ) {

  return value * dimwits::coulombs;
}

inline ChannelRadius toChannelRadius( double value ) {

  return value * njoy::resonanceReconstruction::rootBarns;
}

inline Energy toEnergy( double value ) {

  return value * dimwits::electronVolt;
}

inline ReducedWidth toReducedWidth( double value ) {

  return value * dimwits::rootElectronVolt;
}

inline QValue toQValue( double value ) {

  return value * dimwits::electronVolt;
}

inline std::vector< Energy > toEnergyArray( const std::vector< double > array ) {

  std::vector< Energy > energies;
  energies.reserve( array.size() );
  for ( auto value : array ) {

    energies.push_back( toEnergy( value ) );
  }
  return energies;
}

inline std::vector< ReducedWidth > toReducedWidthArray( const std::vector< double > array ) {

  std::vector< ReducedWidth > widths;
  widths.reserve( array.size() );
  for ( auto value : array ) {

    widths.push_back( toReducedWidth( value ) );
  }
  return widths;
}

#endif
