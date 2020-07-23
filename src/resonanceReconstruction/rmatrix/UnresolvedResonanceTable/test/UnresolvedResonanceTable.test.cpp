#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using UnresolvedResonanceTable = rmatrix::UnresolvedResonanceTable;
using UnresolvedResonance = rmatrix::UnresolvedResonance;
using ChannelID = rmatrix::ChannelID;

#include "resonanceReconstruction/rmatrix/UnresolvedResonanceTable/test/UnresolvedResonanceTable.test.hpp"
#include "resonanceReconstruction/rmatrix/UnresolvedResonanceTable/test/call.test.hpp"
