namespace multiLevel {

using ResonanceShape = decltype( psiChi( 1. * electronVolt ) );

/*
using TemperatureDependentResonanceShape =
  decltype( psiChi( 1. * electronVolt, 1. * kelvin ) );
*/

#include "resonanceReconstruction/breitWigner/multiLevel/src/D.hpp"

#include "resonanceReconstruction/breitWigner/multiLevel/Resonance.hpp"
#include "resonanceReconstruction/breitWigner/multiLevel/Lvalue.hpp"
#include "resonanceReconstruction/breitWigner/multiLevel/Base.hpp"
#include "resonanceReconstruction/breitWigner/multiLevel/Type.hpp"
#include "resonanceReconstruction/breitWigner/multiLevel/Apply.hpp"

}
