add_executable( resonanceReconstruction.rmatrix.legacy.unresolved.SpinGroup.test SpinGroup.test.cpp )
target_link_libraries( resonanceReconstruction.rmatrix.legacy.unresolved.SpinGroup.test PUBLIC resonanceReconstruction )
add_test( NAME resonanceReconstruction.rmatrix.legacy.unresolved.SpinGroup COMMAND resonanceReconstruction.rmatrix.legacy.unresolved.SpinGroup.test )
