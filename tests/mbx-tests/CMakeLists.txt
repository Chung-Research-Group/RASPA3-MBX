add_executable(unit_tests_mbx
  mbx.cpp
  main.cpp)

add_custom_command(TARGET unit_tests_mbx POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          "${CMAKE_CURRENT_SOURCE_DIR}/mbx.json"
          "$<TARGET_FILE_DIR:unit_tests_mbx>/mbx.json"
  COMMENT "Copying mbx.json next to unit_tests_mbx")

raspa_target_link_libraries(unit_tests_mbx)
target_link_libraries(unit_tests_mbx PUBLIC
                      GTest::gtest
                      GTest::gtest_main
                      GTest::gmock
                      GTest::gmock_main)

gtest_add_tests(
  TARGET unit_tests_mbx
  TEST_SUFFIX .noArgs
  TEST_LIST noArgsTests
  WORKING_DIRECTORY "$<TARGET_FILE_DIR:unit_tests_mbx>")
  
set_tests_properties(${noArgsTests}
  PROPERTIES
    ENVIRONMENT "RASPA_MBX_TEST_JSON=$<TARGET_FILE_DIR:unit_tests_mbx>/mbx.json")

if(BUILD_MAC_PACKAGE)
  install(CODE "execute_process(COMMAND
      codesign --force --options runtime --timestamp --sign \"Developer ID Application: David Dubbeldam (24U2ZRZ6SC)\"  
      \"tests/raspakit-tests/unit_tests_mbx\")" COMPONENT unit_tests)
endif()

install(TARGETS unit_tests_mbx DESTINATION share/raspa3/tests COMPONENT unit_tests)