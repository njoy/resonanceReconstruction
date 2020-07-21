#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;
using namespace njoy::resonanceReconstruction::rmatrix;

// convenience typedefs
using OrbitalAngularMomentum = rmatrix::OrbitalAngularMomentum;
using Spin = rmatrix::Spin;
using ResonanceRange = ENDF::ResonanceRange;

constexpr AtomicMass neutronMass = 1.008664 * daltons;
constexpr ElectricalCharge elementaryCharge = 1.602e-19 * coulomb;

#include "resonanceReconstruction/rmatrix/test/possibleChannelTotalAngularMomentumValues.test.hpp"
#include "resonanceReconstruction/rmatrix/test/possibleChannelSpinValues.test.hpp"
#include "resonanceReconstruction/rmatrix/test/calculateCoulombPhaseShift.test.hpp"
#include "resonanceReconstruction/rmatrix/test/fromENDF.test.hpp"
