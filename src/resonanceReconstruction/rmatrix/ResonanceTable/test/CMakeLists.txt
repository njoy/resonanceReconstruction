add_executable( resonanceReconstruction.rmatrix.ResonanceTable.test ResonanceTable.test.cpp )
target_link_libraries( resonanceReconstruction.rmatrix.ResonanceTable.test PUBLIC resonanceReconstruction ) 
add_test( NAME resonanceReconstruction.rmatrix.ResonanceTable COMMAND resonanceReconstruction.rmatrix.ResonanceTable.test )
