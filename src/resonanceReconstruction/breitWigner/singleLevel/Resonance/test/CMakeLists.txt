
add_executable( resonanceReconstruction.breitWigner.singleLevel.Resonance.test Resonance.test.cpp )
target_compile_options( resonanceReconstruction.breitWigner.singleLevel.Resonance.test PRIVATE ${${PREFIX}_common_flags}
$<$<BOOL:${strict}>:${${PREFIX}_strict_flags}>$<$<CONFIG:DEBUG>:
${${PREFIX}_DEBUG_flags}
$<$<BOOL:${coverage}>:${${PREFIX}_coverage_flags}>>
$<$<CONFIG:RELEASE>:
${${PREFIX}_RELEASE_flags}
$<$<BOOL:${link_time_optimization}>:${${PREFIX}_link_time_optimization_flags}>
$<$<BOOL:${nonportable_optimization}>:${${PREFIX}_nonportable_optimization_flags}>>

${CXX_appended_flags} ${resonanceReconstruction_appended_flags} )
target_link_libraries( resonanceReconstruction.breitWigner.singleLevel.Resonance.test PUBLIC resonanceReconstruction ) 
add_test( NAME resonanceReconstruction.breitWigner.singleLevel.Resonance COMMAND resonanceReconstruction.breitWigner.singleLevel.Resonance.test )