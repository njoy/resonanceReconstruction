add_executable( resonanceReconstruction.rmatrix.SpinGroup.test SpinGroup.test.cpp )
target_link_libraries( resonanceReconstruction.rmatrix.SpinGroup.test PUBLIC resonanceReconstruction )
add_test( NAME resonanceReconstruction.rmatrix.SpinGroup COMMAND resonanceReconstruction.rmatrix.SpinGroup.test )
