
add_executable( resonanceReconstruction.breitWigner.singleLevel.Type.test Type.test.cpp )
target_compile_options( resonanceReconstruction.breitWigner.singleLevel.Type.test PRIVATE ${${PREFIX}_common_flags}
$<$<BOOL:${strict}>:${${PREFIX}_strict_flags}>$<$<CONFIG:DEBUG>:
${${PREFIX}_DEBUG_flags}
$<$<BOOL:${coverage}>:${${PREFIX}_coverage_flags}>>
$<$<CONFIG:RELEASE>:
${${PREFIX}_RELEASE_flags}
$<$<BOOL:${link_time_optimization}>:${${PREFIX}_link_time_optimization_flags}>
$<$<BOOL:${nonportable_optimization}>:${${PREFIX}_nonportable_optimization_flags}>>

${CXX_appended_flags} ${resonanceReconstruction_appended_flags} )
target_link_libraries( resonanceReconstruction.breitWigner.singleLevel.Type.test PUBLIC resonanceReconstruction ) 
add_test( NAME resonanceReconstruction.breitWigner.singleLevel.Type COMMAND resonanceReconstruction.breitWigner.singleLevel.Type.test )