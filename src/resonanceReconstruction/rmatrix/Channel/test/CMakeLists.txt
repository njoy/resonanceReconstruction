add_executable( resonanceReconstruction.rmatrix.Channel.test Channel.test.cpp )
target_link_libraries( resonanceReconstruction.rmatrix.Channel.test PUBLIC resonanceReconstruction )
add_test( NAME resonanceReconstruction.rmatrix.Channel COMMAND resonanceReconstruction.rmatrix.Channel.test )
