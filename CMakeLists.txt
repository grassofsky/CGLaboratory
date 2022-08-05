cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(cg_laboratory)

list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(TEXTURE_PATH ${CMAKE_SOURCE_DIR}/resources/textures/)
set(MESHES_PATH ${CMAKE_SOURCE_DIR}/resources/meshes/)
configure_file(
    "${CMAKE_SOURCE_DIR}/include/config.h.in"
    "${CMAKE_SOURCE_DIR}/include/config.h"
)

set(EXTERNAL_INCLUDE ${CMAKE_SOURCE_DIR}/include/external/)

include_directories(
    ${CMAKE_SOURCE_DIR}/include
)

add_subdirectory(src/PPLL_OIT)
add_subdirectory(src/meshdiff)
add_subdirectory(src/LearnOpenGL)