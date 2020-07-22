#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using Resonance = rmatrix::Resonance;

#include "../src/resonanceReconstruction/rmatrix/Resonance/test/Resonance.test.hpp"
#include "../src/resonanceReconstruction/rmatrix/Resonance/test/rmatrix.test.hpp"
