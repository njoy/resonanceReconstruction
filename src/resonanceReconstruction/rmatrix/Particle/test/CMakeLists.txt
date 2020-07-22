add_executable( resonanceReconstruction.rmatrix.Particle.test Particle.test.cpp )
target_link_libraries( resonanceReconstruction.rmatrix.Particle.test PUBLIC resonanceReconstruction ) 
add_test( NAME resonanceReconstruction.rmatrix.Particle COMMAND resonanceReconstruction.rmatrix.Particle.test )
