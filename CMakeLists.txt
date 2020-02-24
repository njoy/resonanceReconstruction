
cmake_minimum_required( VERSION 3.2 )
set( CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Supported configuration types" FORCE )
        
if ( NOT GIT_EXECUTABLE )
    find_package( Git 2.1 )
    if ( NOT GIT_FOUND )
        message( FATAL_ERROR "git installation was not found." )
    endif()
endif()

if( NOT ROOT_DIRECTORY )
    set( ROOT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    if ( NOT fetched_subprojects )
        if ( NOT PYTHON_EXECUTABLE )
            find_package( PythonInterp 3.4 )
            if ( NOT PYTHONINTERP_FOUND )
                message( FATAL_ERROR "Python interpreter installation was not found." )
            endif()
        endif()
        execute_process( COMMAND ${PYTHON_EXECUTABLE} "./metaconfigure/fetch_subprojects.py"
                         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                         RESULT_VARIABLE fetch_failure )
        if ( NOT fetch_failure )
            set( fetched_subprojects TRUE CACHE BOOL "fetch script ran" )
        else()
            message( FATAL_ERROR "Failed to fetch dependencies" )
        endif()
    endif()
endif()

if( NOT is_subproject )
    if( NOT ROOT_DIRECTORY )
        set( ROOT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    endif()
    if ( NOT PYTHON_EXECUTABLE )
        find_package( PythonInterp )
        if ( NOT PYTHONINTERP_FOUND )
            message( FATAL_ERROR "Python interpreter installation was not found." )
        endif()
    endif()
    execute_process( COMMAND ${PYTHON_EXECUTABLE} ./metaconfigure/signature.py ${CMAKE_BINARY_DIR}/signature
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE signature_failure )
    file( READ "${CMAKE_BINARY_DIR}/signature.json" SIGNATURE )
endif()

project( resonanceReconstruction LANGUAGES CXX )
get_directory_property( is_subproject PARENT_DIRECTORY )
include( CMakeDependentOption REQUIRED )

set( resonanceReconstruction_AppleClang_minimum_version 8.0 )
set( resonanceReconstruction_GNU_minimum_version 6.2 )
set( resonanceReconstruction_Clang_minimum_version 3.8 )

if( resonanceReconstruction_${CMAKE_CXX_COMPILER_ID}_minimum_version )
    if( CMAKE_CXX_COMPILER_VERSION AND
        CMAKE_CXX_COMPILER_VERSION VERSION_LESS
        ${resonanceReconstruction_${CMAKE_CXX_COMPILER_ID}_minimum_version} )
        message( FATAL_ERROR "${CMAKE_CXX_COMPILER_ID} version must be greater than ${resonanceReconstruction_${CMAKE_CXX_COMPILER_ID}_minimum_version}" )
    endif()
endif()

# general properties
option( resonanceReconstruction_strict "Compile time warnings are converted to errors" ON )

# binary instrumentation
option( coverage "Enable binary instrumentation to collect test coverage information in the DEBUG configuration" )
option( profile_generate "Enable binary instrumentation to generation execution profiles in the RELEASE configuration which may be used to guide later optimization" )

# additional optimizations
option( link_time_optimization "Enable link time optimization in the RELEASE configuration" )
option( profile_use "In the RELEASE configuration, leverage previously generated execution profile to inform optimization decisions" )
option( nonportable_optimization "Enable optimizations which compromise portability of resulting binary in the RELEASE configuration" )

# libraries and linking
option( static "Statically link component and environment libraries" OFF )
if ( static AND ( "${CMAKE_SYSTEM_NAME}" STREQUAL "Darwin" ) )
    message( FATAL_ERROR "Static binaries not supported on OSX" )
endif()

CMAKE_DEPENDENT_OPTION( static_libraries "Statically link component libraries" OFF "NOT static" ON )

option( unit_tests "Compile the resonanceReconstruction unit tests and integrate with ctest" ON ) 

if ( profile_generate AND profile_use )
    message( FATAL_ERROR "Cannot both generate and use execution profile in the same configuration" )
endif()

if ( profile_generate )
    file( MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/profiling" )
endif()

set( resonanceReconstruction_AppleClang_Darwin_common_flags "-stdlib=libc++" "-Wall" "-Wextra" "-Wpedantic" "-std=c++17" )
set( resonanceReconstruction_AppleClang_Darwin_DEBUG_flags "-O0" "-g" "-gdwarf-3" )
set( resonanceReconstruction_AppleClang_Darwin_RELEASE_flags "-O3" "-DNDEBUG" )
set( resonanceReconstruction_AppleClang_Darwin_strict_flags "-Werror" )
set( resonanceReconstruction_AppleClang_Darwin_coverage_flags "--coverage" )
set( resonanceReconstruction_AppleClang_Darwin_subproject_flags  )
set( resonanceReconstruction_AppleClang_Darwin_base_project_flags  )
set( resonanceReconstruction_AppleClang_Darwin_profile_generate_flags "-fprofile-generate='${CMAKE_BINARY_DIR}/profiling'" )
set( resonanceReconstruction_AppleClang_Darwin_link_time_optimization_flags "-flto" )
set( resonanceReconstruction_AppleClang_Darwin_profile_use_flags "-fprofile-use='${CMAKE_BINARY_DIR}/profiling'" )
set( resonanceReconstruction_AppleClang_Darwin_nonportable_optimization_flags "-march=native" )
set( resonanceReconstruction_AppleClang_Darwin_static_flags "-static" )
set( resonanceReconstruction_GNU_Windows_common_flags "-Wno-subobject-linkage" "-Wall" "-Wextra" "-Wpedantic" "-std=c++17" )
set( resonanceReconstruction_GNU_Windows_DEBUG_flags "-O0" "-g" "-gdwarf-3" "-fsignaling-nans" )
set( resonanceReconstruction_GNU_Windows_RELEASE_flags "-O3" "-DNDEBUG" )
set( resonanceReconstruction_GNU_Windows_strict_flags "-Werror" )
set( resonanceReconstruction_GNU_Windows_coverage_flags "--coverage" )
set( resonanceReconstruction_GNU_Windows_subproject_flags  )
set( resonanceReconstruction_GNU_Windows_base_project_flags  )
set( resonanceReconstruction_GNU_Windows_profile_generate_flags "-fprofile-generate='${CMAKE_BINARY_DIR}/profiling'" )
set( resonanceReconstruction_GNU_Windows_link_time_optimization_flags "-flto" )
set( resonanceReconstruction_GNU_Windows_profile_use_flags "-fprofile-use='${CMAKE_BINARY_DIR}/profiling'" )
set( resonanceReconstruction_GNU_Windows_nonportable_optimization_flags "-march=native" )
set( resonanceReconstruction_GNU_Windows_static_flags "-static" )
set( resonanceReconstruction_GNU_Darwin_common_flags "-Wno-subobject-linkage" "-Wall" "-Wextra" "-Wpedantic" "-std=c++17" )
set( resonanceReconstruction_GNU_Darwin_DEBUG_flags "-O0" "-g" "-gdwarf-3" "-fsignaling-nans" )
set( resonanceReconstruction_GNU_Darwin_RELEASE_flags "-O3" "-DNDEBUG" )
set( resonanceReconstruction_GNU_Darwin_strict_flags "-Werror" )
set( resonanceReconstruction_GNU_Darwin_coverage_flags "--coverage" )
set( resonanceReconstruction_GNU_Darwin_subproject_flags  )
set( resonanceReconstruction_GNU_Darwin_base_project_flags  )
set( resonanceReconstruction_GNU_Darwin_profile_generate_flags "-fprofile-generate='${CMAKE_BINARY_DIR}/profiling'" )
set( resonanceReconstruction_GNU_Darwin_link_time_optimization_flags "-flto" )
set( resonanceReconstruction_GNU_Darwin_profile_use_flags "-fprofile-use='${CMAKE_BINARY_DIR}/profiling'" )
set( resonanceReconstruction_GNU_Darwin_nonportable_optimization_flags "-march=native" )
set( resonanceReconstruction_GNU_Darwin_static_flags "-static" )
set( resonanceReconstruction_GNU_MinGW_common_flags "-Wno-subobject-linkage" "-Wall" "-Wextra" "-Wpedantic" "-std=c++17" )
set( resonanceReconstruction_GNU_MinGW_DEBUG_flags "-O0" "-g" "-gdwarf-3" "-fsignaling-nans" )
set( resonanceReconstruction_GNU_MinGW_RELEASE_flags "-O3" "-DNDEBUG" )
set( resonanceReconstruction_GNU_MinGW_strict_flags "-Werror" )
set( resonanceReconstruction_GNU_MinGW_coverage_flags "--coverage" )
set( resonanceReconstruction_GNU_MinGW_subproject_flags  )
set( resonanceReconstruction_GNU_MinGW_base_project_flags  )
set( resonanceReconstruction_GNU_MinGW_profile_generate_flags "-fprofile-generate='${CMAKE_BINARY_DIR}/profiling'" )
set( resonanceReconstruction_GNU_MinGW_link_time_optimization_flags "-flto" )
set( resonanceReconstruction_GNU_MinGW_profile_use_flags "-fprofile-use='${CMAKE_BINARY_DIR}/profiling'" )
set( resonanceReconstruction_GNU_MinGW_nonportable_optimization_flags "-march=native" )
set( resonanceReconstruction_GNU_MinGW_static_flags "-static" )
set( resonanceReconstruction_GNU_CYGWIN_common_flags "-Wno-subobject-linkage" "-Wall" "-Wextra" "-Wpedantic" "-std=gnu++17" )
set( resonanceReconstruction_GNU_CYGWIN_DEBUG_flags "-O0" "-g" "-gdwarf-3" "-fsignaling-nans" )
set( resonanceReconstruction_GNU_CYGWIN_RELEASE_flags "-O3" "-DNDEBUG" )
set( resonanceReconstruction_GNU_CYGWIN_strict_flags "-Werror" )
set( resonanceReconstruction_GNU_CYGWIN_coverage_flags "--coverage" )
set( resonanceReconstruction_GNU_CYGWIN_subproject_flags  )
set( resonanceReconstruction_GNU_CYGWIN_base_project_flags  )
set( resonanceReconstruction_GNU_CYGWIN_profile_generate_flags "-fprofile-generate='${CMAKE_BINARY_DIR}/profiling'" )
set( resonanceReconstruction_GNU_CYGWIN_link_time_optimization_flags "-flto" )
set( resonanceReconstruction_GNU_CYGWIN_profile_use_flags "-fprofile-use='${CMAKE_BINARY_DIR}/profiling'" )
set( resonanceReconstruction_GNU_CYGWIN_nonportable_optimization_flags "-march=native" )
set( resonanceReconstruction_GNU_CYGWIN_static_flags "-static" )
set( resonanceReconstruction_GNU_Linux_common_flags "-Wno-subobject-linkage" "-Wall" "-Wextra" "-Wpedantic" "-std=c++17" )
set( resonanceReconstruction_GNU_Linux_DEBUG_flags "-O0" "-g" "-gdwarf-3" "-fsignaling-nans" )
set( resonanceReconstruction_GNU_Linux_RELEASE_flags "-O3" "-DNDEBUG" )
set( resonanceReconstruction_GNU_Linux_strict_flags "-Werror" )
set( resonanceReconstruction_GNU_Linux_coverage_flags "--coverage" )
set( resonanceReconstruction_GNU_Linux_subproject_flags  )
set( resonanceReconstruction_GNU_Linux_base_project_flags  )
set( resonanceReconstruction_GNU_Linux_profile_generate_flags "-fprofile-generate='${CMAKE_BINARY_DIR}/profiling'" )
set( resonanceReconstruction_GNU_Linux_link_time_optimization_flags "-flto" )
set( resonanceReconstruction_GNU_Linux_profile_use_flags "-fprofile-use='${CMAKE_BINARY_DIR}/profiling'" )
set( resonanceReconstruction_GNU_Linux_nonportable_optimization_flags "-march=native" )
set( resonanceReconstruction_GNU_Linux_static_flags "-static" )
set( resonanceReconstruction_Clang_Windows_common_flags "-stdlib=libc++" "-Wall" "-Wextra" "-Wpedantic" "-std=c++17" )
set( resonanceReconstruction_Clang_Windows_DEBUG_flags "-O0" "-g" "-gdwarf-3" )
set( resonanceReconstruction_Clang_Windows_RELEASE_flags "-O3" "-DNDEBUG" )
set( resonanceReconstruction_Clang_Windows_strict_flags "-Werror" )
set( resonanceReconstruction_Clang_Windows_coverage_flags "--coverage" )
set( resonanceReconstruction_Clang_Windows_subproject_flags  )
set( resonanceReconstruction_Clang_Windows_base_project_flags  )
set( resonanceReconstruction_Clang_Windows_profile_generate_flags "-fprofile-generate='${CMAKE_BINARY_DIR}/profiling'" )
set( resonanceReconstruction_Clang_Windows_link_time_optimization_flags "-flto" )
set( resonanceReconstruction_Clang_Windows_profile_use_flags "-fprofile-use='${CMAKE_BINARY_DIR}/profiling'" )
set( resonanceReconstruction_Clang_Windows_nonportable_optimization_flags "-march=native" )
set( resonanceReconstruction_Clang_Windows_static_flags "-static" )
set( resonanceReconstruction_Clang_Darwin_common_flags "-stdlib=libc++" "-Wall" "-Wextra" "-Wpedantic" "-std=c++17" )
set( resonanceReconstruction_Clang_Darwin_DEBUG_flags "-O0" "-g" "-gdwarf-3" )
set( resonanceReconstruction_Clang_Darwin_RELEASE_flags "-O3" "-DNDEBUG" )
set( resonanceReconstruction_Clang_Darwin_strict_flags "-Werror" )
set( resonanceReconstruction_Clang_Darwin_coverage_flags "--coverage" )
set( resonanceReconstruction_Clang_Darwin_subproject_flags  )
set( resonanceReconstruction_Clang_Darwin_base_project_flags  )
set( resonanceReconstruction_Clang_Darwin_profile_generate_flags "-fprofile-generate='${CMAKE_BINARY_DIR}/profiling'" )
set( resonanceReconstruction_Clang_Darwin_link_time_optimization_flags "-flto" )
set( resonanceReconstruction_Clang_Darwin_profile_use_flags "-fprofile-use='${CMAKE_BINARY_DIR}/profiling'" )
set( resonanceReconstruction_Clang_Darwin_nonportable_optimization_flags "-march=native" )
set( resonanceReconstruction_Clang_Darwin_static_flags "-static" )
set( resonanceReconstruction_Clang_MinGW_common_flags "-Wno-subobject-linkage" "-Wall" "-Wextra" "-Wpedantic" "-std=c++17" )
set( resonanceReconstruction_Clang_MinGW_DEBUG_flags "-O0" "-g" "-gdwarf-3" "-fsignaling-nans" )
set( resonanceReconstruction_Clang_MinGW_RELEASE_flags "-O3" "-DNDEBUG" )
set( resonanceReconstruction_Clang_MinGW_strict_flags "-Werror" )
set( resonanceReconstruction_Clang_MinGW_coverage_flags "--coverage" )
set( resonanceReconstruction_Clang_MinGW_subproject_flags  )
set( resonanceReconstruction_Clang_MinGW_base_project_flags  )
set( resonanceReconstruction_Clang_MinGW_profile_generate_flags "-fprofile-generate='${CMAKE_BINARY_DIR}/profiling'" )
set( resonanceReconstruction_Clang_MinGW_link_time_optimization_flags "-flto" )
set( resonanceReconstruction_Clang_MinGW_profile_use_flags "-fprofile-use='${CMAKE_BINARY_DIR}/profiling'" )
set( resonanceReconstruction_Clang_MinGW_nonportable_optimization_flags "-march=native" )
set( resonanceReconstruction_Clang_MinGW_static_flags "-static" )
set( resonanceReconstruction_Clang_CYGWIN_common_flags "-Wno-subobject-linkage" "-Wall" "-Wextra" "-Wpedantic" "-std=gnu++17" )
set( resonanceReconstruction_Clang_CYGWIN_DEBUG_flags "-O0" "-g" "-gdwarf-3" "-fsignaling-nans" )
set( resonanceReconstruction_Clang_CYGWIN_RELEASE_flags "-O3" "-DNDEBUG" )
set( resonanceReconstruction_Clang_CYGWIN_strict_flags "-Werror" )
set( resonanceReconstruction_Clang_CYGWIN_coverage_flags "--coverage" )
set( resonanceReconstruction_Clang_CYGWIN_subproject_flags  )
set( resonanceReconstruction_Clang_CYGWIN_base_project_flags  )
set( resonanceReconstruction_Clang_CYGWIN_profile_generate_flags "-fprofile-generate='${CMAKE_BINARY_DIR}/profiling'" )
set( resonanceReconstruction_Clang_CYGWIN_link_time_optimization_flags "-flto" )
set( resonanceReconstruction_Clang_CYGWIN_profile_use_flags "-fprofile-use='${CMAKE_BINARY_DIR}/profiling'" )
set( resonanceReconstruction_Clang_CYGWIN_nonportable_optimization_flags "-march=native" )
set( resonanceReconstruction_Clang_CYGWIN_static_flags "-static" )
set( resonanceReconstruction_Clang_Linux_common_flags "-stdlib=libstdc++" "-Wall" "-Wextra" "-Wpedantic" "-std=c++17" )
set( resonanceReconstruction_Clang_Linux_DEBUG_flags "-O0" "-g" "-gdwarf-3" )
set( resonanceReconstruction_Clang_Linux_RELEASE_flags "-O3" "-DNDEBUG" )
set( resonanceReconstruction_Clang_Linux_strict_flags "-Werror" )
set( resonanceReconstruction_Clang_Linux_coverage_flags "--coverage" )
set( resonanceReconstruction_Clang_Linux_subproject_flags  )
set( resonanceReconstruction_Clang_Linux_base_project_flags  )
set( resonanceReconstruction_Clang_Linux_profile_generate_flags "-fprofile-generate='${CMAKE_BINARY_DIR}/profiling'" )
set( resonanceReconstruction_Clang_Linux_link_time_optimization_flags "-flto" )
set( resonanceReconstruction_Clang_Linux_profile_use_flags "-fprofile-use='${CMAKE_BINARY_DIR}/profiling'" )
set( resonanceReconstruction_Clang_Linux_nonportable_optimization_flags "-march=native" )
set( resonanceReconstruction_Clang_Linux_static_flags "-static" )
get_directory_property( is_subproject PARENT_DIRECTORY )

if( NOT TARGET fmt-adapter )
    add_subdirectory( ${ROOT_DIRECTORY}/subprojects/fmt-adapter )
endif()

if( NOT TARGET range-v3-adapter )
    add_subdirectory( ${ROOT_DIRECTORY}/subprojects/range-v3-adapter )
endif()

if( NOT TARGET spdlog-adapter )
    add_subdirectory( ${ROOT_DIRECTORY}/subprojects/spdlog-adapter )
endif()

if( NOT TARGET hana-adapter )
    add_subdirectory( ${ROOT_DIRECTORY}/subprojects/hana-adapter )
endif()

if( NOT TARGET catch-adapter )
    add_subdirectory( ${ROOT_DIRECTORY}/subprojects/catch-adapter )
endif()

if( NOT TARGET Log )
    add_subdirectory( ${ROOT_DIRECTORY}/subprojects/Log )
endif()

if( NOT TARGET disco )
    add_subdirectory( ${ROOT_DIRECTORY}/subprojects/disco )
endif()

if( NOT TARGET header-utilities )
    add_subdirectory( ${ROOT_DIRECTORY}/subprojects/header-utilities )
endif()

if( NOT TARGET ENDFtk )
    add_subdirectory( ${ROOT_DIRECTORY}/subprojects/ENDFtk )
endif()

if( NOT TARGET eigen-adapter )
    add_subdirectory( ${ROOT_DIRECTORY}/subprojects/eigen-adapter )
endif()

if( NOT TARGET dimwits )
    add_subdirectory( ${ROOT_DIRECTORY}/subprojects/dimwits )
endif()

if( NOT TARGET interpolation )
    add_subdirectory( ${ROOT_DIRECTORY}/subprojects/interpolation )
endif()


if ( NOT GIT_EXECUTABLE )
    find_package( Git )
    if ( NOT GIT_FOUND )
        message( FATAL_ERROR "git installation was not found." )
    endif()
endif()

execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
) 

message( STATUS "" )
message( STATUS "-----------------------------------------------------------" )
message( STATUS "" )
message( STATUS "resonanceReconstruction" )
message( STATUS "Git current branch: ${GIT_BRANCH}" )
message( STATUS "Git commit hash: ${GIT_HASH}" )
message( STATUS "" )
message( STATUS "-----------------------------------------------------------" )

add_library( resonanceReconstruction INTERFACE )
target_sources( resonanceReconstruction INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/EnergyRange.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/pack.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/ZeroWidth.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/CrossSection.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/CRTP.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/lvalue.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/multilevel.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/resonance.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/singleLevel.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/Type.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/CRTP/src/call.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/CRTP/src/derived.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/CRTP/src/evaluate.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/CRTP/src/isTemperature.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/lvalue/Type.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/lvalue/Type/src/call.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/lvalue/Type/src/competitiveWidth.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/lvalue/Type/src/ctor.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/lvalue/Type/src/penetrationShift.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/lvalue/Type/src/phaseShift.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/lvalue/Type/src/statisticalFactorSummation.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/lvalue/Type/src/waveNumber.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/lvalue/Type/src/withCompetitiveWidth.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/lvalue/Type/src/withoutCompetitiveWidth.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/multilevel/Apply.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/multilevel/Base.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/multilevel/Lvalue.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/multilevel/Resonance.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/multilevel/Type.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/multilevel/Apply/src/build.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/multilevel/Apply/src/call.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/multilevel/Apply/src/competitiveWidth.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/multilevel/Apply/src/resonance.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/multilevel/Base/src/ctor.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/multilevel/Base/src/evaluate.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/multilevel/Base/src/lvalues.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/multilevel/Base/src/order.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/multilevel/Base/src/psiChi.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/multilevel/Lvalue/src/call.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/multilevel/Lvalue/src/competitiveWidth.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/multilevel/Lvalue/src/D.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/multilevel/Lvalue/src/evaluate.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/multilevel/Lvalue/src/resonances.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/multilevel/Resonance/src/call.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/multilevel/src/D.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/multilevel/src/psiChi.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/multilevel/Type/src/ctor1.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/multilevel/Type/src/ctor2.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/multilevel/Type/src/evaluate1.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/multilevel/Type/src/evaluate2.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/resonance/Type.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/resonance/src/call.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/resonance/Type/src/ctor.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/singleLevel/Apply.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/singleLevel/Base.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/singleLevel/Lvalue.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/singleLevel/Resonance.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/singleLevel/Type.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/singleLevel/Apply/src/build.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/singleLevel/Apply/src/call.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/singleLevel/Apply/src/competitiveWidth.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/singleLevel/Apply/src/lvalue.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/singleLevel/Apply/src/resonance.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/singleLevel/Base/src/evaluate.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/singleLevel/Base/src/lvalues.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/singleLevel/Base/src/psiChi.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/singleLevel/Lvalue/src/call.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/singleLevel/Lvalue/src/competitiveWidth.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/singleLevel/Lvalue/src/evaluate.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/singleLevel/Lvalue/src/evaluated.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/singleLevel/Lvalue/src/resonances.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/singleLevel/Resonance/src/call.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/singleLevel/Type/src/ctor1.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/singleLevel/Type/src/ctor2.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/singleLevel/Type/src/evaluate1.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/singleLevel/Type/src/evaluate2.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/src/psiChi.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/Type/src/ctor.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/Type/src/isTemperature.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/breitWigner/Type/src/neutronWaveNumber.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Apply.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Base.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Inspector.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Lvalue.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Resonance.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Type.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Apply/src/build.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Apply/src/call.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Apply/src/lvalue.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Apply/src/resonance.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Base/src/ctor.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Base/src/evaluate.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Base/src/neutronWaveNumber.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Base/src/order.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Lvalue/src/ap.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Lvalue/src/call.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Lvalue/src/ctor.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Lvalue/src/hasFission.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Lvalue/src/penetrationShift.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Lvalue/src/phaseShift.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Lvalue/src/solveRfunction.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Lvalue/src/solveRmatrix.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Resonance/src/call.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Resonance/src/ctor.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Type/src/ctor1.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Type/src/ctor2.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Type/src/evaluate1.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/reichMoore/Type/src/evaluate2.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/src/channelRadius.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/src/neutronWaveNumber.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/src/penetrationShift.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/src/phaseShift.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/src/radius.hpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/resonanceReconstruction/src/root.hpp" )

target_include_directories( resonanceReconstruction INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/src src )

set( PREFIX resonanceReconstruction_${CMAKE_CXX_COMPILER_ID}_${CMAKE_SYSTEM_NAME} )

target_link_libraries( resonanceReconstruction INTERFACE "$<$<AND:$<CONFIG:RELEASE>,$<BOOL:${link_time_optimization}>>:${${PREFIX}_RELEASE_flags};${${PREFIX}_link_time_optimization_flags}$<$<BOOL:${profile_generate}>:${${PREFIX}_profile_generate_flags};>$<$<BOOL:${profile_use}>:${${PREFIX}_profile_use_flags};>$<$<BOOL:${nonportable_optimization}>:${${PREFIX}_nonportable_optimization_flags};>>$<$<CONFIG:DEBUG>:$<$<BOOL:${coverage}>:${${PREFIX}_coverage_flags};>>$<$<BOOL:CXX_appended_flags>:${CXX_appended_flags};>$<$<BOOL:resonanceReconstruction_appended_flags>:${resonanceReconstruction_appended_flags};>" )

target_link_libraries( resonanceReconstruction INTERFACE dimwits INTERFACE interpolation INTERFACE ENDFtk INTERFACE eigen-adapter )

if( NOT is_subproject )
    enable_testing() 
    if ( unit_tests )
        add_subdirectory( src/resonanceReconstruction/breitWigner/singleLevel/Resonance/test )
        add_subdirectory( src/resonanceReconstruction/reichMoore/Type/test )
        add_subdirectory( src/resonanceReconstruction/breitWigner/singleLevel/Apply/test )
        add_subdirectory( src/resonanceReconstruction/reichMoore/Apply/test )
        add_subdirectory( src/resonanceReconstruction/breitWigner/multilevel/Type/test )
        add_subdirectory( src/resonanceReconstruction/breitWigner/multilevel/Resonance/test )
        add_subdirectory( src/resonanceReconstruction/breitWigner/singleLevel/Lvalue/test )
        add_subdirectory( src/resonanceReconstruction/breitWigner/singleLevel/Type/test )
        add_subdirectory( src/resonanceReconstruction/breitWigner/multilevel/Apply/test )
        add_subdirectory( src/resonanceReconstruction/test )
    endif() 
endif()

install( FILES "${CMAKE_BINARY_DIR}/signature.json"
         DESTINATION share/${CMAKE_PROJECT_NAME}
         PERMISSIONS OWNER_READ OWNER_WRITE
                     GROUP_READ
                     WORLD_READ
)
install( DIRECTORY src/ DESTINATION include
         FILE_PERMISSIONS OWNER_READ OWNER_WRITE
                          GROUP_READ
                          WORLD_READ
         DIRECTORY_PERMISSIONS OWNER_READ OWNER_WRITE
                               GROUP_READ
                               WORLD_READ
         FILES_MATCHING REGEX ".*\.h[+][+]$|.*\.hxx$|.*\.hpp$|.*\.hh$|.*\.h$"  )

                    INCLUDE(CPack)
                    