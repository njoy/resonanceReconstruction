add_executable( resonanceReconstruction.rmatrix.legacy.resolved.SpinGroup.MultiLevelBreitWigner.test MultiLevelBreitWigner.test.cpp )
target_link_libraries( resonanceReconstruction.rmatrix.legacy.resolved.SpinGroup.MultiLevelBreitWigner.test PUBLIC resonanceReconstruction )
add_test( NAME resonanceReconstruction.rmatrix.legacy.resolved.SpinGroup.MultiLevelBreitWigner COMMAND resonanceReconstruction.rmatrix.legacy.resolved.SpinGroup.MultiLevelBreitWigner.test )
