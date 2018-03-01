cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0022 NEW)

project(histo)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

FILE(COPY histo.hpp DESTINATION include)
install(FILES histo.hpp DESTINATION include)

set(ENABLE_GOOGLE_TEST "OFF" CACHE BOOL "Requires Gtest installed")

IF(ENABLE_GOOGLE_TEST)
    add_subdirectory(test)
ENDIF()
