add_executable( resonanceReconstruction.rmatrix.test rmatrix.test.cpp )
target_link_libraries( resonanceReconstruction.rmatrix.test PUBLIC resonanceReconstruction )
add_test( NAME resonanceReconstruction.rmatrix COMMAND resonanceReconstruction.rmatrix.test )
