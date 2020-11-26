#ifndef NJOY_R2_QUANTITIES
#define NJOY_R2_QUANTITIES

// system includes

// other includes
#include "dimwits.hpp"

namespace njoy {
namespace resonanceReconstruction {

  using namespace dimwits;

  // physical units
  using RootBarn = decltype( unit::sqrt( Barns() ) );
  using RootBarns = RootBarn;
  using InvRootBarn = decltype( pow( RootBarn(), Ratio<-1> ) );
  using InvRootBarns = InvRootBarn;
  using RootElectronVolt = decltype( unit::sqrt( ElectronVolts() ) );
  using RootElectronVolts = RootElectronVolt;
  using InvElectronVolt = decltype( pow( ElectronVolts(), Ratio<-1> ) );
  using InvElectronVolts = InvElectronVolt;
  using ElectronVoltSecond = decltype( ElectronVolts() * Seconds() );
  using ElectronVoltSeconds = ElectronVoltSecond;
  using ElectronVoltSquared = decltype( ElectronVolts() * ElectronVolts() );
  using FaradPerMeter = decltype( farad / meter );
  using CoulombSquaredSecondPerMeter = decltype( Coulomb() * Coulomb() *
                                                 Seconds() / Meters() );

  // constants
  constexpr double pi = 3.141592653589793;

  // constants - taken from 2014 CODATA
  constexpr Quantity< ElectronVoltSecond > hbar = 6.582119514e-16 * electronVolt * second;
  constexpr Quantity< FaradPerMeter > epsilon0 = 8.854187817e-12 * farad / meter;

  // unit variables
  constexpr Quantity< RootBarn > rootBarn = 1.0E-12 * centi(meter);
  constexpr Quantity< RootBarn > rootBarns = rootBarn;
  constexpr Quantity< RootElectronVolt > rootElectronVolt = 1.0 * unit::sqrt( electronVolt );
  constexpr Quantity< RootElectronVolt > rootElectronVolts = rootElectronVolt;

  // physical quantities
  using AtomicMass = Quantity< Dalton >;
  using CrossSection = Quantity< Barn >;
  using ElectricalCharge = Quantity< Coulomb >;
  using Energy = Quantity< ElectronVolt >;
  using EnergySquared = Quantity< ElectronVoltSquared >;
  using LevelSpacing = Quantity< ElectronVolt >;
  using QValue = Quantity< ElectronVolt >;
  using WaveNumber = Quantity< InvRootBarn >;
  using ChannelRadius = Quantity< RootBarn >;
  using EtaParameter = Quantity< CoulombSquaredSecondPerMeter >;
  using Width = Quantity< ElectronVolt >;
  using ReducedWidth = Quantity< RootElectronVolt >;
  using FluctuationIntegral = Quantity< InvElectronVolt >;

} // resonanceReconstruction namespace
} // njoy namespace

#endif
