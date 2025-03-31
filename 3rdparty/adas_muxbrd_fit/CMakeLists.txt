# CMakeList.txt : CMake project for adas_muxbrd_sync include source and define
# project specific logic here.
#
#cmake_minimum_required(VERSION 3.8)
#
#project("adas_muxbrd_fit")
#
#find_package(ROOT REQUIRED)
#include_directories(${PROJECT_SOURCE_DIR}/inc ${PROJECT_SOURCE_DIR}/../../my_pub/lib ${PROJECT_SOURCE_DIR}/../../my_pub/common_lib ${ROOT_INCLUDE_DIR})
#include(${ROOT_USE_FILE})
#
#file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/../../my_pub/src/*.cc)
#file(GLOB headers ${PROJECT_SOURCE_DIR}/inc/*.hh ${PROJECT_SOURCE_DIR}/../../my_pub/lib/*.hh ${PROJECT_SOURCE_DIR}/../../my_pub/common_lib)
#
## Add source to this project's executable.
#add_executable(adas_muxbrd_sync "adas_muxbrd_fit.cc" "adas_muxbrd_fit.hh" ${sources} ${headers})
#target_link_libraries(adas_muxbrd_sync ${ROOT_LIBRARIES})
#target_link_libraries(adas_muxbrd_sync Spectrum)

# TODO: Add tests and install targets if needed.

set(exec_name "adas_muxbrd_fit")
set(aget_code_path ${PROJECT_SOURCE_DIR}/3rdparty)
set(this_path "adas_muxbrd_fit")
set(pub_path "my_pub")
find_package(Eigen3 REQUIRED)

file(GLOB sources
  ${aget_code_path}/${this_path}/adas_muxbrd_fit.cc
  ${aget_code_path}/${this_path}/src/*
  ${aget_code_path}/${pub_path}/src/*
  )

add_executable(${exec_name} ${sources})
target_include_directories(
  ${exec_name}
  PUBLIC
  ${ROOT_INCLUDE_DIRS}
  ${aget_code_path}/${this_path}
  ${aget_code_path}/${this_path}/inc
  ${aget_code_path}/${pub_path}/lib
  ${aget_code_path}/${pub_path}/common_lib
)

file(GLOB config_files
  ${aget_code_path}/${this_path}/calc/*
  )
foreach(x ${config_files})
  file (COPY ${x} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file (COPY ${x} DESTINATION ${CMAKE_BINARY_DIR}/calc)
endforeach()

target_link_libraries(
  ${exec_name}
  ${ROOT_LIBRARIES}
  ROOT::Spectrum
  )
