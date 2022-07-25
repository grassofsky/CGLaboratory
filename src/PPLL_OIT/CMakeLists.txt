project(ppll_oit)

find_package(OpenSceneGraph REQUIRED osg osgViewer osgUtil)

include_directories(
    ${OPENSCENEGRAPH_INCLUDE_DIRS}
)

set(PROJECT_SOURCE_FILE_LIST
    main.cpp
)

add_executable(${PROJECT_NAME} ${PROJECT_SOURCE_FILE_LIST})
target_link_libraries(${PROJECT_NAME} ${OPENSCENEGRAPH_LIBRARIES})
