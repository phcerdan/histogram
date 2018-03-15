find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(test_histo test_histo.cpp)
target_link_libraries(test_histo histo)
target_link_libraries(test_histo ${GTEST_BOTH_LIBRARIES})

if(WITH_VTK)
add_executable(test_visualize_histo test_visualize_histo.cpp)
target_link_libraries(test_visualize_histo histo)
target_link_libraries(test_visualize_histo ${GTEST_BOTH_LIBRARIES})
endif()
