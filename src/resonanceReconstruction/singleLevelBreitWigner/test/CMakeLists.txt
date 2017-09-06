
add_executable( resonanceReconstruction.singleLevelBreitWigner.test
                channelRadius.test.cpp
                neutronWaveNumber.test.cpp
                penetrationShift.test.cpp
                phaseShift.test.cpp
                radius.test.cpp
                singleLevelBreitWigner.test.cpp )
target_compile_options( resonanceReconstruction.singleLevelBreitWigner.test PRIVATE ${${PREFIX}_common_flags}
$<$<BOOL:${strict}>:${${PREFIX}_strict_flags}>$<$<CONFIG:DEBUG>:
${${PREFIX}_DEBUG_flags}
$<$<BOOL:${coverage}>:${${PREFIX}_coverage_flags}>>
$<$<CONFIG:RELEASE>:
${${PREFIX}_RELEASE_flags}
$<$<BOOL:${link_time_optimization}>:${${PREFIX}_link_time_optimization_flags}>
$<$<BOOL:${nonportable_optimization}>:${${PREFIX}_nonportable_optimization_flags}>>

${CXX_appended_flags} ${resonanceReconstruction_appended_flags} )
target_link_libraries( resonanceReconstruction.singleLevelBreitWigner.test PUBLIC resonanceReconstruction ) 
add_test( NAME resonanceReconstruction.singleLevelBreitWigner COMMAND resonanceReconstruction.singleLevelBreitWigner.test )