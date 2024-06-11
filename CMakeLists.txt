########################################################################
# Preamble
########################################################################
cmake_minimum_required( VERSION 3.14 )

set( subproject OFF )
if( DEFINED PROJECT_NAME )
  set( subproject ON )
endif()

project( resonanceReconstruction
  VERSION 0.1.0
  LANGUAGES CXX
)

include( CTest )
include( CMakeDependentOption )
include( GNUInstallDirs )

########################################################################
# Project-wide setup
########################################################################


set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_STANDARD_REQUIRED YES )

cmake_dependent_option(
  resonanceReconstruction.tests
  "Build the resonanceReconstruction unit tests and integrate with ctest" ON
  "BUILD_TESTING AND NOT ${subproject}" OFF
)

########################################################################
# Dependencies
########################################################################

set( REPOSITORIES "release"
    CACHE STRING
    "Options for where to fetch repositories: develop, release, local"
    )

message( STATUS "Using ${REPOSITORIES} repositories" )

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

add_library( resonanceReconstruction INTERFACE )
add_library( njoy::resonanceReconstruction ALIAS resonanceReconstruction )
target_include_directories( resonanceReconstruction
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/src>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries( resonanceReconstruction
  INTERFACE
    ENDFtk
    eigen
    interpolation
    dimwits
    elementary
)

if( CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR )
  if( resonanceReconstruction.tests )
    include( cmake/unit_testing.cmake )
  endif()
endif()
