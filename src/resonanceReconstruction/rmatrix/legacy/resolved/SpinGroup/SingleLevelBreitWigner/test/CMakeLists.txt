add_executable( resonanceReconstruction.rmatrix.legacy.resolved.SpinGroup.SingleLevelBreitWigner.test SingleLevelBreitWigner.test.cpp )
target_link_libraries( resonanceReconstruction.rmatrix.legacy.resolved.SpinGroup.SingleLevelBreitWigner.test PUBLIC resonanceReconstruction )
add_test( NAME resonanceReconstruction.rmatrix.legacy.resolved.SpinGroup.SingleLevelBreitWigner COMMAND resonanceReconstruction.rmatrix.legacy.resolved.SpinGroup.SingleLevelBreitWigner.test )
