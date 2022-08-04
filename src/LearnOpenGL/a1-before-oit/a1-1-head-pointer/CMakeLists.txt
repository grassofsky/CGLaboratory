project(a1-1-head-pointer)

find_package(OpenSceneGraph REQUIRED osg osgViewer osgUtil osgDB)
find_package(OpenMesh)

include_directories(
    ${OPENSCENEGRAPH_INCLUDE_DIRS}
    ${OPENMESH_INCLUDE_DIRS}
)

set(PROJECT_SOURCE_FILE_LIST
    head-pointer.cpp
    head-pointer.vert
    head-pointer.frag
)

add_executable(${PROJECT_NAME} ${PROJECT_SOURCE_FILE_LIST})
target_link_libraries(${PROJECT_NAME} 
    ${OPENSCENEGRAPH_LIBRARIES}
    ${OPENMESH_LIBRARIES})
