add_executable( resonanceReconstruction.rmatrix.LMatrixCalculator.ShiftFactor.test ShiftFactor.test.cpp )
target_link_libraries( resonanceReconstruction.rmatrix.LMatrixCalculator.ShiftFactor.test PUBLIC resonanceReconstruction ) 
add_test( NAME resonanceReconstruction.rmatrix.LMatrixCalculator.ShiftFactor COMMAND resonanceReconstruction.rmatrix.LMatrixCalculator.ShiftFactor.test )
