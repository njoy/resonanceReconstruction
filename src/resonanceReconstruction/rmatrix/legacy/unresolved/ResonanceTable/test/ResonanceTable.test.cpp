// include Catch2
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
using Catch::Approx;

// what we are testing
#include "resonanceReconstruction/rmatrix/legacy/unresolved/ResonanceTable.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using ResonanceTable = rmatrix::legacy::unresolved::ResonanceTable;
using Resonance = rmatrix::legacy::unresolved::Resonance;
using Degrees = rmatrix::legacy::unresolved::Degrees;

#include "resonanceReconstruction/rmatrix/legacy/unresolved/ResonanceTable/test/ResonanceTable.test.hpp"
#include "resonanceReconstruction/rmatrix/legacy/unresolved/ResonanceTable/test/call.test.hpp"
