#ifndef NJOY_RESONANCE_RECONSTRUCTION
#define NJOY_RESONANCE_RECONSTRUCTION

#include <complex>
#include <iterator>
#include <memory>

#include "interpolation.hpp"
#include "dimwits.hpp"
#include "elementary.hpp"

#include "range/v3/distance.hpp"
#include "range/v3/iterator_range.hpp"
#include "range/v3/to_container.hpp"
#include "range/v3/action/sort.hpp"
#include "range/v3/action/unique.hpp"
#include "range/v3/algorithm/count_if.hpp"
#include "range/v3/algorithm/for_each.hpp"
#include "range/v3/algorithm/find_if.hpp"
#include "range/v3/numeric/accumulate.hpp"
#include "range/v3/utility/iterator.hpp"
#include "range/v3/view/all.hpp"
#include "range/v3/view/cartesian_product.hpp"
#include "range/v3/view/chunk.hpp"
#include "range/v3/view/concat.hpp"
#include "range/v3/view/filter.hpp"
#include "range/v3/view/group_by.hpp"
#include "range/v3/view/indices.hpp"
#include "range/v3/view/iota.hpp"
#include "range/v3/view/join.hpp"
#include "range/v3/view/repeat_n.hpp"
#include "range/v3/view/single.hpp"
#include "range/v3/view/transform.hpp"
#include "range/v3/view/zip.hpp"
#include "range/v3/view/zip_with.hpp"

// ENDF components
#include "ENDFtk/section/2/151.hpp"
namespace endf {

  using ScatteringRadius = njoy::ENDFtk::section::Type< 2, 151 >::ScatteringRadius;
  using ResonanceRange = njoy::ENDFtk::section::Type< 2, 151 >::ResonanceRange;
  using BreitWignerLValue = njoy::ENDFtk::section::Type< 2, 151 >::BreitWignerLValue;
  using ReichMooreLValue = njoy::ENDFtk::section::Type< 2, 151 >::ReichMooreLValue;
  using SingleLevelBreitWigner = njoy::ENDFtk::section::Type< 2, 151 >::SingleLevelBreitWigner;
  using MultiLevelBreitWigner = njoy::ENDFtk::section::Type< 2, 151 >::MultiLevelBreitWigner;
  using ReichMoore = njoy::ENDFtk::section::Type< 2, 151 >::ReichMoore;
  using RMatrixLimited = njoy::ENDFtk::section::Type< 2, 151 >::RMatrixLimited;
  using UnresolvedEnergyIndependent = njoy::ENDFtk::section::Type< 2, 151 >::UnresolvedEnergyIndependent;
  using UnresolvedEnergyDependentFissionWidths = njoy::ENDFtk::section::Type< 2, 151 >::UnresolvedEnergyDependentFissionWidths;
  using UnresolvedEnergyDependent = njoy::ENDFtk::section::Type< 2, 151 >::UnresolvedEnergyDependent;
}

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
using CrossSection = Quantity< Barn >;
using ElectricalCharge = Quantity< Coulomb >;
using Energy = Quantity< ElectronVolt >;
using LevelSpacing = Quantity< ElectronVolt >;
using QValue = Quantity< ElectronVolt >;
using WaveNumber = Quantity< InvRootBarn >;
using ChannelRadius = Quantity< RootBarn >;
using EtaParameter = Quantity< CoulombSquaredSecondPerMeter >;
using Width = Quantity< ElectronVolt >;
using ReducedWidth = Quantity< RootElectronVolt >;
using FluctuationIntegral = Quantity< InvElectronVolt >;

#include "resonanceReconstruction/rmatrix.hpp"
}
}

#endif
