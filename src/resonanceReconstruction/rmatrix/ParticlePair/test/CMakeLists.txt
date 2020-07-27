add_executable( resonanceReconstruction.rmatrix.ParticlePair.test ParticlePair.test.cpp )
target_link_libraries( resonanceReconstruction.rmatrix.ParticlePair.test PUBLIC resonanceReconstruction )
add_test( NAME resonanceReconstruction.rmatrix.ParticlePair COMMAND resonanceReconstruction.rmatrix.ParticlePair.test )
