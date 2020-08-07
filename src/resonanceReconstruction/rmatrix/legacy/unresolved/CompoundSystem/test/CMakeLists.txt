add_executable( resonanceReconstruction.rmatrix.legacy.unresolved.CompoundSystem.test CompoundSystem.test.cpp )
target_link_libraries( resonanceReconstruction.rmatrix.legacy.unresolved.CompoundSystem.test PUBLIC resonanceReconstruction )
add_test( NAME resonanceReconstruction.rmatrix.legacy.unresolved.CompoundSystem COMMAND resonanceReconstruction.rmatrix.legacy.unresolved.CompoundSystem.test )
