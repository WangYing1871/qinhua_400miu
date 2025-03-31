#cmake_minimum_required (VERSION 3.8)
#
#project ("calcu_and_draw")
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
#	${PROJECT_SOURCE_DIR}/../../my_pub/common_lib/*.hh
#)
#
#
## Add source to this project's executable.
#
#add_executable (calcu_and_draw "calcu_and_draw.cc" "calcu_and_draw.hh" ${sources} ${headers})
#target_link_libraries(calcu_and_draw ${ROOT_LIBRARIES})

set(exec_name "calcu_and_draw")
set(aget_code_path ${PROJECT_SOURCE_DIR}/3rdparty)
set(this_path "calcu_and_draw")
set(pub_path "my_pub")
file(GLOB sources
  ${aget_code_path}/${this_path}/calcu_and_draw.cc
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
target_link_libraries(
  ${exec_name}
  ${ROOT_LIBRARIES}
  ROOT::Spectrum
  )
#file(GLOB config_files
#  ${aget_code_path}/${this_path}/calc/*
#  )
#foreach(x ${config_files})
#  file (COPY ${x} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
#  file (COPY ${x} DESTINATION ${CMAKE_BINARY_DIR}/calc)
#endforeach()
set(install_path ${CMAKE_BINARY_DIR}/bin/aget)
install(TARGETS ${exec_name} RUNTIME DESTINATION ${install_path})
