#######################################################################
# Setup
#######################################################################

message( STATUS "Adding resonanceReconstruction unit testing" )
enable_testing()


#######################################################################
# Unit testing directories
#######################################################################

add_subdirectory( src/resonanceReconstruction/rmatrix/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/Particle/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/ParticlePair/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/ChannelQuantumNumbers/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/ChannelRadii/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/Channel/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/BaseResonance/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/Resonance/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/UnresolvedResonance/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/ResonanceTable/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/LMatrixCalculator/ShiftFactor/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/LMatrixCalculator/Constant/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/SpinGroup/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/CompoundSystem/test )

add_subdirectory( src/resonanceReconstruction/rmatrix/legacy/unresolved/test )
add_subdirectory( src/resonanceReconstruction/rmatrix/legacy/unresolved/Data/test )

add_subdirectory( src/resonanceReconstruction/test )
add_subdirectory( src/resonanceReconstruction/reichMoore/Apply/test )
add_subdirectory( src/resonanceReconstruction/reichMoore/Type/test )
add_subdirectory( src/resonanceReconstruction/breitWigner/multiLevel/Apply/test )
add_subdirectory( src/resonanceReconstruction/breitWigner/multiLevel/Type/test )
add_subdirectory( src/resonanceReconstruction/breitWigner/multiLevel/Resonance/test )
add_subdirectory( src/resonanceReconstruction/breitWigner/singleLevel/Apply/test )
add_subdirectory( src/resonanceReconstruction/breitWigner/singleLevel/Lvalue/test )
add_subdirectory( src/resonanceReconstruction/breitWigner/singleLevel/Type/test )
add_subdirectory( src/resonanceReconstruction/breitWigner/singleLevel/Resonance/test )
