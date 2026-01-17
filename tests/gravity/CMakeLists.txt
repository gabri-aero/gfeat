file(GLOB TEST_FILES "Test*.cpp")


foreach(TEST_FILE ${TEST_FILES})
    get_filename_component(TEST_FILE_NAME ${TEST_FILE} NAME_WE)
    add_executable(${TEST_FILE_NAME} ${TEST_FILE})
    target_link_libraries(${TEST_FILE_NAME} PRIVATE ${CPP_TARGET} GTest::gtest_main)
    target_include_directories(${TEST_FILE_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/external/fft)
    target_include_directories(${TEST_FILE_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/external/functions)
    target_include_directories(${TEST_FILE_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
    set_target_properties(${TEST_FILE_NAME} PROPERTIES
    INSTALL_RPATH "$ORIGIN"
    BUILD_WITH_INSTALL_RPATH TRUE
    )
    gtest_discover_tests(${TEST_FILE_NAME})
    message(STATUS "${TEST_FILE_NAME} added as an executable")
endforeach()