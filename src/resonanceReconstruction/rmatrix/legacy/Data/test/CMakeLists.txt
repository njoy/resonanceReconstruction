add_executable( resonanceReconstruction.rmatrix.legacy.unresolved.Data.test Data.test.cpp )
target_link_libraries( resonanceReconstruction.rmatrix.legacy.unresolved.Data.test PUBLIC resonanceReconstruction )
add_test( NAME resonanceReconstruction.rmatrix.legacy.unresolved.Data COMMAND resonanceReconstruction.rmatrix.legacy.unresolved.Data.test )
