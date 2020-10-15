add_executable( resonanceReconstruction.rmatrix.legacy.Data.test Data.test.cpp )
target_link_libraries( resonanceReconstruction.rmatrix.legacy.Data.test PUBLIC resonanceReconstruction )
add_test( NAME resonanceReconstruction.rmatrix.legacy.Data COMMAND resonanceReconstruction.rmatrix.legacy.Data.test )
