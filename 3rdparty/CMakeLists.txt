add_subdirectory(fec2det_dec)
add_subdirectory(hit_position_reconstrcut)
add_subdirectory(tracker_seek)
add_subdirectory(adas_muxbrd_fit)
#add_subdirectory(calcu_and_draw)
#option(BUILD_CALC_AND_DRAW "Build calc_and_draw" OFF)
#option(BUILD_CALC_AND_DRAW "Build calc_and_draw" ON)
if (BUILD_CALC_AND_DRAW)
  message(STATUS "build calcu_and_draw enable")
  add_subdirectory(calcu_and_draw)
endif()
