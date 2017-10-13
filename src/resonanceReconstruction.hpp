#ifndef NJOY_RESONANCE_RECONSTRUCTION
#define NJOY_RESONANCE_RECONSTRUCTION

#include <complex>

#include "ENDFtk.hpp"
#include "interpolation.hpp"
#include "dimwits.hpp"

#include <Eigen/Dense>

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

}
}

#endif
