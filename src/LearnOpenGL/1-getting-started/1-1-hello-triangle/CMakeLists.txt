project(1_1_hello_triangle)

find_package(OpenSceneGraph REQUIRED osg osgViewer osgUtil osgDB)

include_directories(
    ${OPENSCENEGRAPH_INCLUDE_DIRS}
)

set(PROJECT_SOURCE_FILE_LIST
    hello_triangle.cpp
)

add_executable(${PROJECT_NAME} ${PROJECT_SOURCE_FILE_LIST})
target_link_libraries(${PROJECT_NAME} ${OPENSCENEGRAPH_LIBRARIES})
