find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
set(tests_)

add_executable(test_histo test_histo.cpp)
target_link_libraries(test_histo histo)
target_link_libraries(test_histo ${GTEST_BOTH_LIBRARIES})
list(APPEND tests_ test_histo)

if(WITH_VTK)
add_executable(test_visualize_histo test_visualize_histo.cpp)
target_link_libraries(test_visualize_histo histo)
target_link_libraries(test_visualize_histo ${GTEST_BOTH_LIBRARIES})
list(APPEND tests_ test_visualize_histo)
endif()

gtest_discover_tests(
    ${tests_}
    TEST_PREFIX ${SG_MODULE_NAME}||${test_name}||
    PROPERTIES LABELS ${SG_MODULE_NAME}
    )
