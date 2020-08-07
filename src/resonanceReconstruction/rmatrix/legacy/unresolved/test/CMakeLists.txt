add_executable( resonanceReconstruction.rmatrix.legacy.unresolved.test unresolved.test.cpp )
target_link_libraries( resonanceReconstruction.rmatrix.legacy.unresolved.test PUBLIC resonanceReconstruction )
add_test( NAME resonanceReconstruction.rmatrix.legacy.unresolved COMMAND resonanceReconstruction.rmatrix.legacy.unresolved.test )
