add_executable( resonanceReconstruction.rmatrix.legacy.resolved.SpinGroup.test SpinGroup.test.cpp )
target_link_libraries( resonanceReconstruction.rmatrix.legacy.resolved.SpinGroup.test PUBLIC resonanceReconstruction )
add_test( NAME resonanceReconstruction.rmatrix.legacy.resolved.SpinGroup COMMAND resonanceReconstruction.rmatrix.legacy.resolved.SpinGroup.test )
