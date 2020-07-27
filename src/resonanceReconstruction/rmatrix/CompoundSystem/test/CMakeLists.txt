add_executable( resonanceReconstruction.rmatrix.CompoundSystem.test CompoundSystem.test.cpp )
target_link_libraries( resonanceReconstruction.rmatrix.CompoundSystem.test PUBLIC resonanceReconstruction )
add_test( NAME resonanceReconstruction.rmatrix.CompoundSystem COMMAND resonanceReconstruction.rmatrix.CompoundSystem.test )
