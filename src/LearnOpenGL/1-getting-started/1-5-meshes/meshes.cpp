/// see: https://learnopengl-cn.github.io/01%20Getting%20started/06%20Textures/
///
/// ��Ҫ���ƵĶ����ܹ��ڴ����Ͻ�����ת���ţ���Ҫ��shader��ʹ��osg_ModelViewProjectionMatrix
/// ͬʱ    viewer.getCamera()->getGraphicsContext()->getState()->setUseModelViewAndProjectionUniforms(true);
/// ��������osg_ModelViewProjectionMatrix
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
