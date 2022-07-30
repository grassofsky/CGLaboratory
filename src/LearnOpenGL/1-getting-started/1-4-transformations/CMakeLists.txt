project(1-4-transformations)

find_package(OpenSceneGraph REQUIRED osg osgViewer osgUtil osgDB osgGA)

include_directories(
    ${EXTERNAL_INCLUDE}
    ${OPENSCENEGRAPH_INCLUDE_DIRS}
)

set(PROJECT_SOURCE_FILE_LIST
    transformations.cpp
)

add_executable(${PROJECT_NAME} ${PROJECT_SOURCE_FILE_LIST})
target_link_libraries(${PROJECT_NAME} ${OPENSCENEGRAPH_LIBRARIES})
