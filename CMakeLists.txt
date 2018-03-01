cmake_minimum_required(VERSION 2.8)
set(CMAKE_ECLIPSE_VERSION=4.4)
cmake_policy(SET CMP0022 NEW)

project(histo-header)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

add_library(histo-header SHARED histo-header.h)
set_target_properties(histo-header PROPERTIES LINKER_LANGUAGE CXX)
FILE(COPY histo-header.h DESTINATION include)

# TEST
### Google TEST AND MOCK
# There must be an environment variable (in .bashrc or profile.d)
# export GMOCK_HOME=/path/to/gmock-1.7.0
set(ENABLE_GOOGLE_TEST "OFF" CACHE BOOL "Requires ENV variable GMOCK_HOME pointing to gmock parent folder")

IF(ENABLE_GOOGLE_TEST)
    add_subdirectory(test)
ENDIF()
