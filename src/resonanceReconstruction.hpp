#ifndef NJOY_RESONANCE_RECONSTRUCTION
#define NJOY_RESONANCE_RECONSTRUCTION

#include "ENDFtk.hpp"
#include "interpolation.hpp"
#include "dimwits.hpp"

namespace njoy {
namespace resonanceReconstruction {

using namespace dimwits;

using RootBarn = decltype( unit::sqrt( Barns() ) );
using RootBarns = RootBarn;
using InvRootBarn = decltype( pow( RootBarn(), Ratio<-1> ) );
using InvRootBarns = InvRootBarn;
using InvElectronVolt = decltype( pow( ElectronVolts(), Ratio<-1> ) );
using InvElectronVolts = InvElectronVolt;

constexpr Quantity< RootBarn > rootBarn = 1.0E-12 * centi(meter);
constexpr  Quantity< RootBarn > rootBarns = rootBarn;

constexpr double pi = 3.141592653589793;

#include "resonanceReconstruction/singleLevelBreitWigner.hpp" 

}
}

#endif
