#ifndef NJOY_RESONANCE_RECONSTRUCTION
#define NJOY_RESONANCE_RECONSTRUCTION

#include <complex>

#include "ENDFtk.hpp"
#include "interpolation.hpp"
#include "dimwits.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wextra"
#include <Eigen/Dense>
#pragma GCC diagnostic pop

namespace njoy {
namespace resonanceReconstruction {

using namespace dimwits;

constexpr double pi = 3.141592653589793;

using Matrix3x3 = Eigen::Matrix3cd;
using RootBarn = decltype( unit::sqrt( Barns() ) );
using RootBarns = RootBarn;
using InvRootBarn = decltype( pow( RootBarn(), Ratio<-1> ) );
using InvRootBarns = InvRootBarn;
using RootElectronVolt = decltype( unit::sqrt( ElectronVolts() ) );
using RootElectronVolts = RootElectronVolt;
using InvElectronVolt = decltype( pow( ElectronVolts(), Ratio<-1> ) );
using InvElectronVolts = InvElectronVolt;

constexpr Quantity< RootBarn > rootBarn = 1.0E-12 * centi(meter);
constexpr Quantity< RootBarn > rootBarns = rootBarn;

template< int i >
using Integer = std::integral_constant< int, i >;

namespace ENDF = ENDFtk::resonanceParameters;

#include "resonanceReconstruction/src/radius.hpp"
#include "resonanceReconstruction/src/channelRadius.hpp"
#include "resonanceReconstruction/src/root.hpp"
#include "resonanceReconstruction/src/neutronWaveNumber.hpp"
#include "resonanceReconstruction/src/penetrationShift.hpp"
#include "resonanceReconstruction/src/phaseShift.hpp"

#include "resonanceReconstruction/pack.hpp"
#include "resonanceReconstruction/EnergyRange.hpp"
#include "resonanceReconstruction/ZeroWidth.hpp"

#include "resonanceReconstruction/breitWigner.hpp"
#include "resonanceReconstruction/reichMoore.hpp"

// physical constants and quantities

// physical units
using ElectronVoltSecond = decltype( ElectronVolts() * Seconds() );
using ElectronVoltSeconds = ElectronVoltSecond;
using CoulombSquaredSecondPerMeter = decltype( Coulomb() * Coulomb() *
                                               Seconds() / Meters() );
using FaradPerMeter = decltype( farad / meter );

// hbar constant in eV s - taken from 2014 CODATA
constexpr Quantity< ElectronVoltSecond > hbar = 6.582119514e-16 * electronVolt * second;

// epsilon0 constant in F m^-1 - taken from 2014 CODATA
constexpr Quantity< FaradPerMeter > epsilon0 = 8.854187817e-12 * farad / meter;

constexpr Quantity< RootElectronVolt > rootElectronVolt = 1.0 * unit::sqrt( electronVolt );
constexpr Quantity< RootElectronVolt > rootElectronVolts = rootElectronVolt;

// physical quantities
using AtomicMass = Quantity< Dalton >;
using ElectricalCharge = Quantity< Coulomb >;
using Energy = Quantity< ElectronVolt >;
using QValue = Quantity< ElectronVolt >;
using WaveNumber = Quantity< InvRootBarn >;
using ChannelRadius = Quantity< RootBarn >;
using EtaParameter = Quantity< CoulombSquaredSecondPerMeter >;
using ReducedWidth = Quantity< RootElectronVolt >;

#include "resonanceReconstruction/rmatrix.hpp"
}
}

#endif
