add_executable( resonanceReconstruction.rmatrix.Resonance.test Resonance.test.cpp )
target_link_libraries( resonanceReconstruction.rmatrix.Resonance.test PUBLIC resonanceReconstruction )
add_test( NAME resonanceReconstruction.rmatrix.Resonance COMMAND resonanceReconstruction.rmatrix.Resonance.test )
