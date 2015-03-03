set(GMOCK_HOME "$ENV{GMOCK_HOME}" CACHE PATH "Set this to gmock parent folder")
SET( Boost_LIBRARIES ${Boost_LIBRARIES} PARENT_SCOPE )
#FIND GMOCK
find_path(GTEST_INCLUDE_DIR
         gtest/gtest.h
         HINTS ${GMOCK_HOME}/gtest/include
         PATH_SUFFIXES gtest
         DOC "gtest include dir")
find_path(GMOCK_INCLUDE_DIR
         gmock/gmock.h
         HINTS ${GMOCK_HOME}/include
         PATH_SUFFIXES gmock
         DOC "gmock include dir")
include_directories(${GMOCK_INCLUDE_DIR} ${GTEST_INCLUDE_DIR})

find_library(GTEST_LIBRARY
            NAMES gtest
                  libgtest
                  libgtest.a
            HINTS ${GMOCK_HOME}/mybuild/gtest
                  ${GMOCK_HOME}/gtest
                  ${GMOCK_HOME}/gtest/mybuild)
find_library(GTEST_LIBRARY_MAIN
            NAMES gtest_main
                  libgtest_main
                  libgtest.a
            HINTS ${GMOCK_HOME}/mybuild/gtest
                  ${GMOCK_HOME}/gtest
                  ${GMOCK_HOME}/gtest/mybuild)

find_library(GMOCK_LIBRARY
            NAMES gmock
                  libgmock
                  libgmock.a
            HINTS ${GMOCK_HOME}/mybuild
                  ${GMOCK_HOME})
find_library(GMOCK_LIBRARY_MAIN
            NAMES gmock_main
                  libgmock_main
                  libgmock.a
            HINTS ${GMOCK_HOME}/mybuild
                  ${GMOCK_HOME})

include_directories("..")
find_package(Threads)
set_source_files_properties (test_histo.cpp PROPERTIES COMPILE_FLAGS "-Wno-sign-compare")
add_executable(test_histo.exe test_histo.cpp)
target_link_libraries(test_histo.exe ${GMOCK_LIBRARY_MAIN} ${GTEST_LIBRARY_MAIN})
target_link_libraries(test_histo.exe ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(test_histo.exe histo-header)

if(WITH_BOOST)
add_executable(test_histo-boost.exe test_histo-boost.cpp )
target_link_libraries(test_histo-boost.exe ${GMOCK_LIBRARY_MAIN} ${GTEST_LIBRARY_MAIN})
target_link_libraries(test_histo-boost.exe ${CMAKE_THREAD_LIBS_INIT})
# target_link_libraries(test_histo-boost.exe ${PROJECT_BINARY_DIR}/lib/libhisto-header-boost.a)
target_link_libraries(test_histo-boost.exe  ${Boost_Libraries})
# add_library(test-library-boost ../histo-header-boost.h)
# target_link_libraries(test-library-boost PUBLIC ${Boost_Libraries})
# set_target_properties(test-library-boost PROPERTIES LINKER_LANGUAGE CXX)
# target_link_libraries(test_histo-boost.exe test-library-boost)
# target_link_libraries(test_histo-boost.exe histo-header-boost)
# SET( Boost_Libraries_DIR ${Boost_Libraries_DIR} PARENT_SCOPE)
# message("with boost?  ${WITH_BOOST}")
# IF(WITH_BOOST)
#     SET( Boost_Libraries ${Boost_Libraries} PARENT_SCOPE)
#     target_link_libraries(test_histo.exe ${Boost_Libraries})
# endif()
endif()
