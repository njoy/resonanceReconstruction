#######################################################################
# Setup
#######################################################################

message( STATUS "Adding resonanceReconstruction unit testing" )
enable_testing()

FetchContent_MakeAvailable( Catch2 )

function( add_cpp_test name source )

  set( test_name "resonanceReconstruction.${name}.test" )
  add_executable( ${test_name} ${source} )
  add_test( NAME ${test_name} COMMAND ${test_name} )
  target_link_libraries( ${test_name} PRIVATE resonanceReconstruction )
  target_link_libraries( ${test_name} PRIVATE Catch2::Catch2WithMain )

  file( GLOB resources "resources/*" )
  foreach( resource ${resources} )
    file( COPY "${resource}" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}" )
  endforeach()

endfunction()

#######################################################################
# Unit testing directories
#######################################################################

add_subdirectory( src/resonanceReconstruction/breitWigner/multiLevel/Apply/test )
add_subdirectory( src/resonanceReconstruction/breitWigner/multiLevel/Resonance/test )
add_subdirectory( src/resonanceReconstruction/breitWigner/multiLevel/Type/test )
add_subdirectory( src/resonanceReconstruction/breitWigner/singleLevel/Apply/test )
add_subdirectory( src/resonanceReconstruction/breitWigner/singleLevel/Lvalue/test )
add_subdirectory( src/resonanceReconstruction/breitWigner/singleLevel/Resonance/test )
add_subdirectory( src/resonanceReconstruction/breitWigner/singleLevel/Type/test )
add_subdirectory( src/resonanceReconstruction/reichMoore/Apply/test )
add_subdirectory( src/resonanceReconstruction/reichMoore/Type/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/Channel/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/ChannelQuantumNumbers/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/ChannelRadii/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/ChannelRadiusTable/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/CompoundSystem/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/LMatrixCalculator/Constant/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/LMatrixCalculator/ShiftFactor/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/Particle/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/ParticleChannelData/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/ParticlePair/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/Resonance/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/ResonanceTable/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/RLMatrixCalculator/ReichMoore/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/SpinGroup/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/legacy/Data/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/legacy/resolved/CompoundSystem/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/legacy/resolved/Resonance/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/legacy/resolved/SpinGroup/MultiLevelBreitWigner/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/legacy/resolved/SpinGroup/SingleLevelBreitWigner/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/legacy/unresolved/CompoundSystem/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/legacy/unresolved/Resonance/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/legacy/unresolved/ResonanceTable/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/legacy/unresolved/SpinGroup/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/legacy/unresolved/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/test )
add_subdirectory( src/resonanceReconstruction/test )
