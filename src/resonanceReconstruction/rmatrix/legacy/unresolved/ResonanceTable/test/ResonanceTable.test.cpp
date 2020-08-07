#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using ResonanceTable = rmatrix::legacy::unresolved::ResonanceTable;
using Resonance = rmatrix::legacy::unresolved::Resonance;
using Degrees = rmatrix::legacy::unresolved::Degrees;

#include "resonanceReconstruction/rmatrix/legacy/unresolved/ResonanceTable/test/ResonanceTable.test.hpp"
#include "resonanceReconstruction/rmatrix/legacy/unresolved/ResonanceTable/test/call.test.hpp"
