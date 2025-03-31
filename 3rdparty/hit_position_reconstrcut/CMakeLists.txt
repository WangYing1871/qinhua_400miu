# CMakeList.txt : CMake project for hit_position_reconstrcut, include source and define
# project specific logic here.
#
#cmake_minimum_required(VERSION 3.8)
#
#project("hit_position_reconstrcut")
#
#find_package(ROOT REQUIRED)
#include_directories(
#	${PROJECT_SOURCE_DIR}/inc
#	${PROJECT_SOURCE_DIR}/../../my_pub/lib
#	${PROJECT_SOURCE_DIR}/../../my_pub/common_lib
#	${ROOT_INCLUDE_DIR}
#)
#include(${ROOT_USE_FILE})
#
#file(GLOB sources
#	${PROJECT_SOURCE_DIR}/src/*.cc
#	${PROJECT_SOURCE_DIR}/../../my_pub/src/*.cc
#)
#file(GLOB headers
#	${PROJECT_SOURCE_DIR}/inc/*.hh
#	${PROJECT_SOURCE_DIR}/../../my_pub/lib/*.hh
#	${PROJECT_SOURCE_DIR}/../../my_pub/common_lib
#)
#
## Add source to this project's executable.
#add_executable(hit_position_reconstrcut "hit_position_reconstrcut.cc" "hit_position_reconstrcut.hh" ${sources} ${headers})
#target_link_libraries(hit_position_reconstrcut ${ROOT_LIBRARIES})
#
##set(REQUIRED_COPY_FILES
##	z_posiiton_default.csv
##	link_position_default.csv
##)
#
#foreach(_f ${REQUIRED_COPY_FILES})
#	configure_file(${PROJECT_SOURCE_DIR}/${_f} ${PROJECT_BINARY_DIR}/${_f} COPYONLY)
#endforeach()

set(exec_name "hit_position_reconstrcut")
set(aget_code_path ${PROJECT_SOURCE_DIR}/3rdparty)
set(this_path "hit_position_reconstrcut")
set(pub_path "my_pub")
file(GLOB sources
  ${aget_code_path}/${this_path}/hit_position_reconstrcut.cc
  ${aget_code_path}/${this_path}/src/*
  ${aget_code_path}/${pub_path}/src/*
  )

add_executable(${exec_name} ${sources})
target_include_directories(
  ${exec_name}
  PRIVATE 
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
#message(STATUS "${aget_code_path}/${pub_path}/lib")
#message(STATUS "${aget_code_path}/${pub_path}/common_lib")

target_link_libraries(
  ${exec_name}
  ${ROOT_LIBRARIES}
  )
