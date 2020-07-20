########################################################################
# Preamble
########################################################################

cmake_minimum_required( VERSION 3.14 )
project( resonanceReconstruction LANGUAGES CXX )


########################################################################
# Project-wide setup
########################################################################

set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_STANDARD_REQUIRED YES )

option( resonanceReconstruction_unit_tests
    "Compile the resonanceReconstruction unit tests and integrate with ctest" ON
    )


########################################################################
# Dependencies
########################################################################

include( cmake/resonanceReconstruction_dependencies.cmake )


########################################################################
# Project targets
########################################################################

add_library( resonanceReconstruction INTERFACE )
target_include_directories( resonanceReconstruction INTERFACE src/ )
target_link_libraries( resonanceReconstruction
    INTERFACE ENDFtk
    INTERFACE eigen-adapter
    INTERFACE interpolation
    INTERFACE dimwits
    INTERFACE elementary
    )


#######################################################################
# Top-level Only
#######################################################################

if( CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR )

    # unit testing
    if( resonanceReconstruction_unit_tests )
        include( cmake/unit_testing.cmake )
    endif()

endif()
