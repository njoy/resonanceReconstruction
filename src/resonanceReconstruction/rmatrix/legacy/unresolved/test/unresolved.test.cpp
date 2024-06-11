// include Catch2
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
using Catch::Approx;

// what we are testing
#include "resonanceReconstruction/rmatrix/legacy/unresolved.hpp"

using namespace njoy::resonanceReconstruction;
using namespace njoy::resonanceReconstruction::rmatrix::legacy::unresolved;

// convenience typedefs
using Degrees = rmatrix::legacy::unresolved::Degrees;
using Widths = rmatrix::legacy::unresolved::Widths;
using FluctuationIntegrals = rmatrix::legacy::unresolved::FluctuationIntegrals;

#include "resonanceReconstruction/rmatrix/legacy/unresolved/test/calculateFluctuationIntegrals.test.hpp"
