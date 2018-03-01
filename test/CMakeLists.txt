find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

include_directories("..")
add_executable(test_histo test_histo.cpp)
target_link_libraries(test_histo ${GTEST_BOTH_LIBRARIES})
