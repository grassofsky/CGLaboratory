project(remove_self_inter)

find_package(OpenMesh)

include_directories(
    ${OPENSCENEGRAPH_INCLUDE_DIRS}
    ${OPENMESH_INCLUDE_DIRS}
    ${EXTERNAL_INCLUDE}/GTE
)

set(PROJECT_SOURCE_FILE_LIST
    remove_self_inter.cpp
    box_inter_based_get_self_intersection_pairs.h
    aabsptree_based_get_self_intersection_pairs.h
    ${CGlAB_INCLUDE}/openmesh_ext/openmesh_axis_aligned_bsptree.h
)

add_executable(${PROJECT_NAME} ${PROJECT_SOURCE_FILE_LIST})
target_link_libraries(${PROJECT_NAME} 
    ${OPENMESH_LIBRARIES})
