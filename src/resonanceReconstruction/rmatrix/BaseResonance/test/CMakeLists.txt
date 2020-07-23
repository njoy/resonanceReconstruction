add_executable( resonanceReconstruction.rmatrix.BaseResonance.test BaseResonance.test.cpp )
target_link_libraries( resonanceReconstruction.rmatrix.BaseResonance.test PUBLIC resonanceReconstruction )
add_test( NAME resonanceReconstruction.rmatrix.BaseResonance COMMAND resonanceReconstruction.rmatrix.BaseResonance.test )
