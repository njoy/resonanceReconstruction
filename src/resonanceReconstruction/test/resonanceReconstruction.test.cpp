#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

#include "range/v3/view/linear_distribute.hpp"

using namespace dimwits;
using namespace njoy::resonanceReconstruction;

#include "resonanceReconstruction/test/channelRadius.test.hpp"
#include "resonanceReconstruction/test/radius.test.hpp"
#include "resonanceReconstruction/test/phaseShift.test.hpp"
#include "resonanceReconstruction/test/penetrationShift.test.hpp"
#include "resonanceReconstruction/test/neutronWaveNumber.test.hpp"
