add_executable(unit_tests_mbx energy.cpp)
raspa_target_link_libraries(unit_tests_mbx)
target_link_libraries(unit_tests_mbx PRIVATE raspakit GTest::gtest_main)
target_compile_definitions(unit_tests_mbx PRIVATE
  MBX_SETTINGS_FILE="${CMAKE_SOURCE_DIR}/examples/advanced/8_mbx_co2_in_mfi/mbx.json")
gtest_add_tests(TARGET unit_tests_mbx)

if(BUILD_APP)
  add_test(NAME MBX.RestartRequiresMatchingModel
    COMMAND ${CMAKE_COMMAND} -DRASPA=$<TARGET_FILE:raspa3>
            -DFIXTURE=${CMAKE_SOURCE_DIR}/examples/advanced/8_mbx_co2_in_mfi
            -DWORK=${CMAKE_CURRENT_BINARY_DIR}/restart-model -P ${CMAKE_CURRENT_SOURCE_DIR}/restart.cmake)
endif()
