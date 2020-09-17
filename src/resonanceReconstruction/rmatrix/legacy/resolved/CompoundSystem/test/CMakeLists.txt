add_executable( resonanceReconstruction.rmatrix.legacy.resolved.CompoundSystem.test CompoundSystem.test.cpp )
target_link_libraries( resonanceReconstruction.rmatrix.legacy.resolved.CompoundSystem.test PUBLIC resonanceReconstruction )
add_test( NAME resonanceReconstruction.rmatrix.legacy.resolved.CompoundSystem COMMAND resonanceReconstruction.rmatrix.legacy.resolved.CompoundSystem.test )
