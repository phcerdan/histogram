cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0022 NEW)

project(histo)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(HISTO_HEADERS
    ${INCLUDE_DIR}/histo.hpp
    ${INCLUDE_DIR}/visualize_histo.hpp
    )
# Interface library for header only.
add_library(histo INTERFACE)
target_include_directories(histo INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>  # <prefix>/include/mylib
    )
file(COPY ${HISTO_HEADERS} DESTINATION include)
install(FILES ${HISTO_HEADERS} DESTINATION include)

option(WITH_VTK "Use vtk to visualize histogram" OFF)
if(WITH_VTK)
    find_package(VTK REQUIRED COMPONENTS
        vtkChartsCore
        vtkCommonCore
        vtkCommonDataModel
        vtkInteractionStyle
        vtkRenderingOpenGL2 # avoid vtkRenderer, vtkRenderWindow not found at run time
        vtkRenderingContextOpenGL2 # avoid vtkContextDevice2D not found at run time
        vtkViewsContext2D
        CONFIG)
    include(${VTK_USE_FILE})
    target_link_libraries(histo INTERFACE ${VTK_LIBRARIES})
endif()


set(ENABLE_GOOGLE_TEST "OFF" CACHE BOOL "Requires Gtest installed")

if(ENABLE_GOOGLE_TEST)
    add_subdirectory(test)
endif()
