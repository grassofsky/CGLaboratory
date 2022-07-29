project(1-2-shaders-uniform)

find_package(OpenSceneGraph REQUIRED osg osgViewer osgUtil osgDB)

include_directories(
    ${OPENSCENEGRAPH_INCLUDE_DIRS}
)

set(PROJECT_SOURCE_FILE_LIST
    shader_uniform.cpp
)

add_executable(${PROJECT_NAME} ${PROJECT_SOURCE_FILE_LIST})
target_link_libraries(${PROJECT_NAME} ${OPENSCENEGRAPH_LIBRARIES})
