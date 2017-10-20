namespace multilevel {

using ResonanceShape = decltype( psiChi( 1. * electronVolt ) );

/*
using TemperatureDependentResonanceShape =
  decltype( psiChi( 1. * electronVolt, 1. * kelvin ) );
*/

#include "resonanceReconstruction/breitWigner/multilevel/src/D.hpp"

#include "resonanceReconstruction/breitWigner/multilevel/Resonance.hpp"
#include "resonanceReconstruction/breitWigner/multilevel/Lvalue.hpp"
#include "resonanceReconstruction/breitWigner/multilevel/Base.hpp"
#include "resonanceReconstruction/breitWigner/multilevel/Type.hpp"
#include "resonanceReconstruction/breitWigner/multilevel/Apply.hpp"

}
