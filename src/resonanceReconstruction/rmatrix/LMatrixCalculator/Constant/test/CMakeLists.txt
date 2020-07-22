add_executable( resonanceReconstruction.rmatrix.LMatrixCalculator.Constant.test Constant.test.cpp )
target_link_libraries( resonanceReconstruction.rmatrix.LMatrixCalculator.Constant.test PUBLIC resonanceReconstruction ) 
add_test( NAME resonanceReconstruction.rmatrix.LMatrixCalculator.Constant COMMAND resonanceReconstruction.rmatrix.LMatrixCalculator.Constant.test )
