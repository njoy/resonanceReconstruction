#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;
using namespace njoy::resonanceReconstruction::rmatrix::legacy::unresolved;

// convenience typedefs
using Degrees = rmatrix::legacy::unresolved::Degrees;
using Widths = rmatrix::legacy::unresolved::Widths;
using FluctuationIntegrals = rmatrix::legacy::unresolved::FluctuationIntegrals;

#include "resonanceReconstruction/rmatrix/legacy/unresolved/test/calculateFluctuationIntegrals.test.hpp"
