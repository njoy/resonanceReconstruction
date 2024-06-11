// include Catch2
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
using Catch::Approx;

// what we are testing
#include "resonanceReconstruction/rmatrix/Resonance.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using Resonance = rmatrix::Resonance;

#include "../src/resonanceReconstruction/rmatrix/Resonance/test/Resonance.test.hpp"
