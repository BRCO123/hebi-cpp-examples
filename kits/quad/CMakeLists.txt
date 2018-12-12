cmake_minimum_required (VERSION 3.1)
project (quadruped)

# modified from daisy CMakeList

SET (CMAKE_CXX_STANDARD 11)
SET (CMAKE_CXX_STANDARD_REQUIRED ON)
SET (CMAKE_CXX_FLAGS "-O3")

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

find_package(Qt5 REQUIRED COMPONENTS Core Gui Widgets NO_CMAKE_PATH)
#find_package(Qt5 REQUIRED COMPONENTS Core Gui Widgets 3DCore 3DExtras 3DRender 3DInput NO_CMAKE_PATH)
#find_package(Doxygen)

# Get the C++ API directory:
set(ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../..)
get_filename_component(ROOT_DIR "${ROOT_DIR}" REALPATH)
# `HEBI_CPP_DIR` can be used to point to an already downloaded C++ API.
# This is useful for environments that would prefer not to download files
# from the internet (or even don't have internet access).
if(NOT HEBI_CPP_DIR)
  set(HEBI_DIR ${ROOT_DIR}/hebi-cpp)
  include(../../projects/cmake/DownloadHebiCpp.cmake)
else()
  set(HEBI_DIR ${HEBI_CPP_DIR})
endif()

set(hebi_cpp_build_dir hebi-cpp-${CMAKE_BUILD_TYPE})
string(TOLOWER ${hebi_cpp_build_dir} hebi_cpp_build_dir)
get_filename_component(hebi_cpp_build_dir "${CMAKE_CURRENT_BINARY_DIR}/${hebi_cpp_build_dir}" REALPATH)
add_subdirectory(${HEBI_DIR} ${hebi_cpp_build_dir})

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}/src          # Our source files
  ${HEBI_DIR}/src ${HEBI_DIR}/hebi/include ${HEBI_DIR}/Eigen)
link_directories (
  ${HEBI_CPP_LINK_DIRECTORIES})

# start the different part from Hebi hexapod
add_library(quadruped_core OBJECT
  ${CMAKE_CURRENT_SOURCE_DIR}/src/robot/quadruped_parameters.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/robot/quadruped.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/robot/quadruped_leg.cpp
)

SET(SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/quadruped_control.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/input/input_manager_mobile_io.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/xml_util/pugixml.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/xml_util/xml_helpers.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/util/MadgwickAHRS.cpp
)

add_executable(quadruped_control ${SOURCES} $<TARGET_OBJECTS:quadruped_core>)
qt5_use_modules(quadruped_control Core Gui Widgets)

# link hebi library
if (WIN32) 
target_link_libraries( quadruped_control hebi kernel32 ) # kernel32 for sleep commands.
else()
target_link_libraries( quadruped_control hebi hebic++ m pthread)
endif()

# very fancy tools
# Add ultra-conservative warnings.
if ( CMAKE_COMPILER_IS_GNUCC )
  set_property( TARGET quadruped_control APPEND_STRING PROPERTY COMPILE_FLAGS " -Wall -Wno-int-in-bool-context " )
endif ( CMAKE_COMPILER_IS_GNUCC )
# For Windows, we copy the .dll file into the binary directory so that we don't have to set the PATH variable.
if (WIN32)
add_custom_command( TARGET quadruped_control POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
  "${PROJECT_SOURCE_DIR}/lib/hebi.dll"
  $<TARGET_FILE_DIR:quadruped_control>)
endif()

add_custom_command(TARGET quadruped_control POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                       ${CMAKE_CURRENT_SOURCE_DIR}/resources/ $<TARGET_FILE_DIR:quadruped_control>/)