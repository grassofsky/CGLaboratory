project(meshdiff)

find_package(OpenSceneGraph REQUIRED osg osgViewer osgUtil osgDB osgText)
find_package(OpenMesh)

include_directories(
    ${OPENSCENEGRAPH_INCLUDE_DIRS}
    ${OPENMESH_INCLUDE_DIRS}
    ${EXTERNAL_INCLUDE}/GTE
)

set(PROJECT_SOURCE_FILE_LIST
    meshdiff.cpp
)

add_executable(${PROJECT_NAME} ${PROJECT_SOURCE_FILE_LIST})
target_link_libraries(${PROJECT_NAME} 
    ${OPENSCENEGRAPH_LIBRARIES} 
    ${OPENMESH_LIBRARIES})
