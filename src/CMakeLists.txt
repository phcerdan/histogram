cmake_minimum_required(VERSION 2.8)
set(CMAKE_ECLIPSE_VERSION=4.4)
cmake_policy(SET CMP0022 NEW)
# cmake ../src  -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_ECLIPSE_GENERATE_SOURCE_PROJECT=TRUE
set(WITH_BOOST "OFF" CACHE BOOL "Use boost:filesystem to provide Save() member")
if (UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0  -Wall -std=c++11")
else ()
endif()

project(header-only-project)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)



add_library(histo-header SHARED histo-header.h)
set_target_properties(histo-header PROPERTIES LINKER_LANGUAGE CXX)
FILE(COPY histo-header.h DESTINATION include)

if(WITH_BOOST)
    find_package(Boost COMPONENTS
            system filesystem
            REQUIRED)
    include_directories(${Boost_INCLUDE_DIRS})
    add_library(histo-header-boost SHARED histo-header-boost.h )
    target_link_libraries(histo-header-boost PUBLIC ${Boost_LIBRARIES})
    set_target_properties(histo-header-boost PROPERTIES LINKER_LANGUAGE CXX)
    FILE(COPY histo-header-boost.h DESTINATION include)
endif()

# TEST
### Google TEST AND MOCK
# There must be an environment variable (in .bashrc or profile.d)
# export GMOCK_HOME=/path/to/gmock-1.7.0
set(ENABLE_GOOGLE_TEST "OFF" CACHE BOOL "Requires ENV variable GMOCK_HOME pointing to gmock parent folder")

IF(ENABLE_GOOGLE_TEST)
    add_subdirectory(test)
ENDIF()
