cmake_minimum_required(VERSION 3.1...3.5)
project("400_miu_rad_data")
option (UNIX_USE_CCACHE "Treat warnings as error" ON)
set(CMAKE_CXX_STANDARD 20 CACHE STRING "The C++ standard to use")
set(build_type "Release")
if (CMAKE_BUILD_TYPE)
  set(build_type "${CMAKE_BUILD_TYPE}")
  if ("${build_type}" STREQUAL "Debug")
    add_definitions(-DDEBUG)
  endif()
endif()
add_compile_options(-rdynamic -fPIC)
# CERN-ROOT suggest CXX17
add_compile_options(-Wno-cpp)


message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
find_package(Boost COMPONENTS program_options REQUIRED)
message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")

find_package(ROOT REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/inc
    ${ROOT_INCLUDE_DIRS}
    ${ROOT_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
    )
find_program(ROOTCINT rootcint PATH_SUFFIXES bin REQUIRED)
set(root_cint_out "root_cint_out.cpp")
set(data_struct "data_strcut_cint.h")
execute_process(COMMAND "${ROOTCINT}" -I ${PROJECT_SOURCE_DIR}/inc -f ${root_cint_out} -c ${data_struct} LinkDef.h)
add_subdirectory(src)
if (NOT EXISTS "calc")
  file(COPY "${PROJECT_SOURCE_DIR}/calc" DESTINATION "./")
endif()
add_subdirectory(3rdparty)
