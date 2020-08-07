add_executable( resonanceReconstruction.rmatrix.legacy.unresolved.Resonance.test Resonance.test.cpp )
target_link_libraries( resonanceReconstruction.rmatrix.legacy.unresolved.Resonance.test PUBLIC resonanceReconstruction )
add_test( NAME resonanceReconstruction.rmatrix.legacy.unresolved.Resonance COMMAND resonanceReconstruction.rmatrix.legacy.unresolved.Resonance.test )
