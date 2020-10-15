add_executable( resonanceReconstruction.rmatrix.legacy.resolved.Resonance.test Resonance.test.cpp )
target_link_libraries( resonanceReconstruction.rmatrix.legacy.resolved.Resonance.test PUBLIC resonanceReconstruction )
add_test( NAME resonanceReconstruction.rmatrix.legacy.resolved.Resonance COMMAND resonanceReconstruction.rmatrix.legacy.resolved.Resonance.test )
