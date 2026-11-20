add_executable(unit_tests_mbx energy.cpp)
raspa_target_link_libraries(unit_tests_mbx)
target_link_libraries(unit_tests_mbx PRIVATE raspakit GTest::gtest_main)
target_compile_definitions(unit_tests_mbx PRIVATE
  MBX_SETTINGS_FILE="${CMAKE_SOURCE_DIR}/examples/advanced/8_mbx_co2_in_mfi/mbx.json")
gtest_add_tests(TARGET unit_tests_mbx)
