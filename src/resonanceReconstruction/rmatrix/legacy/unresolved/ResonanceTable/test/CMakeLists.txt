add_executable( resonanceReconstruction.rmatrix.legacy.unresolved.ResonanceTable.test ResonanceTable.test.cpp )
target_link_libraries( resonanceReconstruction.rmatrix.legacy.unresolved.ResonanceTable.test PUBLIC resonanceReconstruction )
add_test( NAME resonanceReconstruction.rmatrix.legacy.unresolved.ResonanceTable COMMAND resonanceReconstruction.rmatrix.legacy.unresolved.ResonanceTable.test )
