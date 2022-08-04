#include <osg/Program>
#include <osg/Geometry>
#include <osg/Geode>
#include <osg/Group>
#include <osgViewer/Viewer>

#include "head-pointer.vert";
#include "head-pointer.frag";

#include "utils/openmesh_to_osg_mesh.h"

int main()
{
    osg::ref_ptr<osg::Vec3Array> vertices = new osg::Vec3Array(3);
    (*vertices)[0].set(-0.5, -0.5, 0.0);
    (*vertices)[1].set(0.5, -0.5, 0.0);
    (*vertices)[2].set(0.0, 0.5, 0.0);

    osg::ref_ptr<osg::DrawElementsUInt> indices = new osg::DrawElementsUInt(GL_TRIANGLES, 3);
    (*indices)[0] = 0;
    (*indices)[1] = 1;
    (*indices)[2] = 2;

    osg::ref_ptr<osg::Vec4Array> colors = new osg::Vec4Array(1);
    (*colors)[0].set(1.0, 0.0, 0.0, 1.0);

    osg::ref_ptr<osg::Geometry> geoTriangle = new osg::Geometry;
    geoTriangle->setDataVariance(osg::Object::DYNAMIC);
    geoTriangle->setUseDisplayList(false);
    geoTriangle->setUseVertexBufferObjects(true);
    geoTriangle->setVertexArray(vertices.get());
    geoTriangle->addPrimitiveSet(indices.get());
    geoTriangle->setColorArray(colors.get());
    geoTriangle->setColorBinding(osg::Geometry::BIND_OVERALL);

    osg::ref_ptr<osg::Geode> geode = new osg::Geode;
    geode->addDrawable(geoTriangle);

    osg::ref_ptr<osg::Program> program = new osg::Program();
    program->addShader(new osg::Shader(osg::Shader::VERTEX, vertexShaderSource));
    program->addShader(new osg::Shader(osg::Shader::FRAGMENT, fragmentShaderSource));

    geode->getOrCreateStateSet()->setAttributeAndModes(program, osg::StateAttribute::ON);

    osg::ref_ptr<osg::Group> root = new osg::Group;
    root->addChild(geode);

    osgViewer::Viewer viewer;
    viewer.setSceneData(root);
    viewer.setUpViewInWindow(100, 100, 800, 600);
    return viewer.run();
}
