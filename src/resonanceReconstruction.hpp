#ifndef NJOY_RESONANCE_RECONSTRUCTION
#define NJOY_RESONANCE_RECONSTRUCTION

#include "ENDFtk.hpp"
#include "interpolation.hpp"
#include "dimwits.hpp"

namespace njoy {
namespace resonanceReconstruction {

constexpr double pi = 3.141592653589793;

using namespace dimwits;

using RootBarn = decltype( unit::sqrt( Barns() ) );
using RootBarns = RootBarn;
using InvRootBarn = decltype( pow( RootBarn(), Ratio<-1> ) );
using InvRootBarns = InvRootBarn;
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

}
}

#endif
