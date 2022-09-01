project(ppll_oit)

find_package(OpenSceneGraph REQUIRED osg osgViewer osgUtil osgDB)
find_package(OpenMesh)

include_directories(
    ${OPENSCENEGRAPH_INCLUDE_DIRS}
    ${OPENMESH_INCLUDE_DIRS}
)

set(PROJECT_SOURCE_FILE_LIST
    build_list_subpass.h
    build_list_subpass.cpp
    resolve_list_subpass.h
    resolve_list_subpass.cpp
    ppll_render_pass.h
    ppll_render_pass.cpp
    ppll_oit.cpp
)

add_executable(${PROJECT_NAME} ${PROJECT_SOURCE_FILE_LIST})
target_link_libraries(${PROJECT_NAME} 
    ${OPENSCENEGRAPH_LIBRARIES} 
    ${OPENMESH_LIBRARIES})
