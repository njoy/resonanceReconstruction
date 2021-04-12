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
option( strict_compile
    "Treat all warnings as errors." ON
    )

# Compile flags
set( common_flags "-Wall" "-Wextra" "-Wpedantic" )
set( strict_flags "-Werror" )
set( release_flags "-O3" )
set( debug_flags "-O0" "-g" )


########################################################################
# Dependencies
########################################################################

set( REPOSITORIES "release"
    CACHE STRING
    "Options for where to fetch repositories: develop, release, local"
    )

if( REPOSITORIES STREQUAL "develop" )
    include( cmake/develop_dependencies.cmake )

elseif( REPOSITORIES STREQUAL "release" )
    include( cmake/release_dependencies.cmake )

elseif( REPOSITORIES STREQUAL "local" )
    include( cmake/local_dependencies.cmake )

endif()


########################################################################
# Project targets
########################################################################

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# resonanceReconstruction : library
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

add_library( resonanceReconstruction INTERFACE
    )
target_include_directories( resonanceReconstruction INTERFACE src/ )
target_link_libraries( resonanceReconstruction
    INTERFACE ENDFtk
    INTERFACE eigen
    INTERFACE interpolation
    INTERFACE dimwits
    INTERFACE elementary
    )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# resonanceReconstruction : python bindings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

include_directories( python/src/ )

pybind11_add_module( resonanceReconstruction.python
    EXCLUDE_FROM_ALL
    python/src/resonanceReconstruction.python.cpp
    python/src/rmatrix/ChannelQuantumNumbers.python.cpp
    python/src/rmatrix/Particle.python.cpp
    )
target_link_libraries( resonanceReconstruction.python PUBLIC resonanceReconstruction )
target_compile_options( resonanceReconstruction.python PUBLIC "-fPIC" )
target_compile_options( resonanceReconstruction.python PUBLIC "-fvisibility=hidden" )
set_target_properties( resonanceReconstruction.python PROPERTIES OUTPUT_NAME resonanceReconstruction )
set_target_properties( resonanceReconstruction.python PROPERTIES COMPILE_DEFINITIONS "PYBIND11" )

#######################################################################
# Top-level Only
#######################################################################

if( CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR )

    # unit testing
    if( resonanceReconstruction_unit_tests )
        include( cmake/unit_testing.cmake )
    endif()

endif()
