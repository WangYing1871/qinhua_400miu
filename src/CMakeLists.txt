function(add_exec sname ename)
  file(GLOB sources
      ${CMAKE_BINARY_DIR}/${root_cint_out}
      ${PROJECT_SOURCE_DIR}/src/${sname}
      ${PROJECT_SOURCE_DIR}/src/util.cpp
      )
  add_executable(${ename} ${sources})
  target_link_libraries(${ename} ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
endfunction()

function(get_file_name path ret)
  string(LENGTH "${path}" length)
  string(FIND "${path}" "/" del_index REVERSE)
  math(EXPR sub_result "${length}-${del_index}")
  math(EXPR start "${del_index}+1")
  string(SUBSTRING "${path}" "${start}" "${sub_result}" name)
  set(${ret} "${name}" PARENT_SCOPE)
endfunction()

include_directories(
    ${PROJECT_SOURCE_DIR}/inc
    ${Boost_INCLUDE_DIRS}
    ${ROOT_INCLUDE_DIRS}

    )
file(GLOB scripts
    ${PROJECT_SOURCE_DIR}/src/*.cpp
    ${PROJECT_SOURCE_DIR}/src/*.cxx
    ${PROJECT_SOURCE_DIR}/src/*.cc
    ${PROJECT_SOURCE_DIR}/src/*.C
    )
set(exec_names PARENT_SCOPE)
foreach(x ${scripts})
  if (NOT "${x}" STREQUAL "${PROJECT_SOURCE_DIR}/src/util.cpp")
    cmake_language(CALL get_file_name "${x}" source_name)
    string(REGEX REPLACE "(\.cpp)|(\.cxx)|(\.cc)|(\.C)" "" exec_name ${source_name})
    cmake_language(CALL add_exec "${source_name}" "${exec_name}")
    list(APPEND exec_names ${exec_name})
  endif()
endforeach()
string(REPLACE ".cpp" "_rdict.pcm" root_cint_out_rdict ${root_cint_out})
if (NOT EXISTS ${CMAKE_BINARY_DIR}/${root_cint_out_rdict})
  cmake_warning(DEPRECATION "${root_cint_out_rdict} can't access! scripts may not exec successful")
else()
  file(COPY
      ${CMAKE_BINARY_DIR}/${root_cint_out_rdict}
      DESTINATION
      "${CMAKE_CURRENT_BINARY_DIR}")
endif()

set(install_des ${CMAKE_BINARY_DIR}/bin)
foreach(x ${exec_names})
  install(TARGETS "${x}"
      RUNTIME DESTINATION ${install_des}
  )
endforeach()
install(FILES ${CMAKE_BINARY_DIR}/${root_cint_out_rdict} DESTINATION ${install_des})
