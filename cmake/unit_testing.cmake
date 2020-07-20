#######################################################################
# Setup
#######################################################################

message( STATUS "Adding NJOY21 unit testing" )
enable_testing()


#######################################################################
# Unit testing directories
#######################################################################

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
