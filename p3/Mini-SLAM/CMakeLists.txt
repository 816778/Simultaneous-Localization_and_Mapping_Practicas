#Miinimum CMAKE version required
cmake_minimum_required(VERSION 3.5)
project(Mini-SLAM)

#Build type:
#   -Set to RELEASE to benchmark your code
#   -Set to DEBUG to debug your code
set(CMAKE_BUILD_TYPE DEBUG)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

#We require the C++14 standard
set(CMAKE_CXX_STANDARD 14)

#Compilation flags:
#   -Flag for code optimizations: -O3
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/Thirdparty/OpenCV/build")
  set(OpenCV_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Thirdparty/OpenCV/build")
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/Thirdparty/Pangolin/build")
  set(Pangolin_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Thirdparty/Pangolin/build")
endif()

add_subdirectory(Lab0)
add_subdirectory(Modules)
add_subdirectory(Apps)
