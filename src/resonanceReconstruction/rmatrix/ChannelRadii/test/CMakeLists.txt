add_executable( resonanceReconstruction.rmatrix.ChannelRadii.test ChannelRadii.test.cpp )
target_link_libraries( resonanceReconstruction.rmatrix.ChannelRadii.test PUBLIC resonanceReconstruction ) 
add_test( NAME resonanceReconstruction.rmatrix.ChannelRadii COMMAND resonanceReconstruction.rmatrix.ChannelRadii.test )
