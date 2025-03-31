#cmake_minimum_required (VERSION 3.8)
#
#project ("fec2det_dec")
#
## set(CMAKE_CXX_STANDARD 11)
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
#add_executable(fec2det_dec "fec2det_dec.cc" "fec2det_dec.hh" ${sources} ${headers})
#target_link_libraries(fec2det_dec ${ROOT_LIBRARIES})
set(exec_name "fec2det_dec")
set(aget_code_path ${PROJECT_SOURCE_DIR}/3rdparty)
set(this_path "fec2det_dec")
set(pub_path "my_pub")
file(GLOB sources
  ${aget_code_path}/${this_path}/fec2det_dec.cc
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

#message(STATUS "${aget_code_path}/${pub_path}/lib")
#message(STATUS "${aget_code_path}/${pub_path}/common_lib")

target_link_libraries(
  ${exec_name}
  ${ROOT_LIBRARIES}
  )
file(GLOB config_files
  ${aget_code_path}/${this_path}/calc/*
  )
foreach(x ${config_files})
  file (COPY ${x} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file (COPY ${x} DESTINATION ${CMAKE_BINARY_DIR}/calc)
endforeach()
set(install_path ${CMAKE_BINARY_DIR}/bin/aget)
install(TARGETS ${exec_name} RUNTIME DESTINATION ${install_path})
