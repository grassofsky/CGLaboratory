/// see: https://learnopengl-cn.github.io/01%20Getting%20started/06%20Textures/
///
/// 想要绘制的对象能够在窗口上进行旋转缩放，需要在shader中使用osg_ModelViewProjectionMatrix
/// 同时    viewer.getCamera()->getGraphicsContext()->getState()->setUseModelViewAndProjectionUniforms(true);
/// 设置启用osg_ModelViewProjectionMatrix
///
#include "config.h"

#include <sstream>

#include <osg/Program>
#include <osg/Geometry>
#include <osg/Geode>
#include <osg/Group>
#include <osgDB/ReadFile>
#include <osgViewer/Viewer>


int main()
{
    osg::ref_ptr<osg::Node> trimesh = osgDB::readNodeFile(std::string(MESHES_PATH) + "sphere.stl");

    osg::ref_ptr<osg::Group> root = new osg::Group;
    root->addChild(trimesh);

    // transform by camera
    osgViewer::Viewer viewer;
    viewer.setSceneData(root);
    viewer.setUpViewInWindow(100, 100, 800, 600);

    return viewer.run();
}
