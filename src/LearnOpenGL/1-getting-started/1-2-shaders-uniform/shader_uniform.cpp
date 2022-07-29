/// see: https://learnopengl-cn.github.io/01%20Getting%20started/05%20Shaders/

#include <osg/Program>
#include <osg/Geometry>
#include <osg/Geode>
#include <osg/Group>
#include <osgViewer/Viewer>

// location = 0 是位置
// location = 3 是颜色
// TODO: 这个关系具体是在哪里体现的
const char* vertexShaderSource = "#version 330 core\n"
"layout (location = 0) in vec4 aPos;\n"
"layout (location = 3) in vec4 aColor;\n"
"out vec4 v2fColor;\n"
"void main()\n"
"{\n"
"   gl_Position = aPos;\n"
"   v2fColor = aColor;\n"
"}\0";
const char* fragmentShaderSource = "#version 330 core\n"
"out vec4 FragColor;\n"
"in vec4 v2fColor;\n"
"uniform vec4 ourColorFactor;\n"
"void main()\n"
"{\n"
"   FragColor = v2fColor * ourColorFactor;\n"
"}\n\0";

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

    osg::ref_ptr<osg::Vec4Array> colors = new osg::Vec4Array(3);
    (*colors)[0].set(0.0, 1.0, 0.0, 1.0);
    (*colors)[1].set(0.0, 1.0, 0.0, 1.0);
    (*colors)[2].set(0.0, 1.0, 0.0, 1.0);

    osg::ref_ptr<osg::Geometry> geoTriangle = new osg::Geometry;
    geoTriangle->setDataVariance(osg::Object::DYNAMIC);
    geoTriangle->setUseDisplayList(false);
    geoTriangle->setUseVertexBufferObjects(true);
    geoTriangle->setUseVertexArrayObject(true);
    geoTriangle->setVertexArray(vertices.get());
    geoTriangle->addPrimitiveSet(indices.get());
    geoTriangle->setColorArray(colors.get());
    geoTriangle->setColorBinding(osg::Geometry::BIND_PER_VERTEX);

    osg::ref_ptr<osg::Geode> geode = new osg::Geode;
    geode->addDrawable(geoTriangle);

    osg::ref_ptr<osg::Program> program = new osg::Program();
    program->addShader(new osg::Shader(osg::Shader::VERTEX, vertexShaderSource));
    program->addShader(new osg::Shader(osg::Shader::FRAGMENT, fragmentShaderSource));

    // set program
    geode->getOrCreateStateSet()->setAttributeAndModes(program, osg::StateAttribute::ON);
    // Set uniform value
    geode->getOrCreateStateSet()->addUniform(new osg::Uniform("ourColorFactor", osg::Vec4f(1.0, 0.5, 0.5, 1.0)));

    osg::ref_ptr<osg::Group> root = new osg::Group;
    root->addChild(geode);

    osgViewer::Viewer viewer;
    viewer.setSceneData(root);
    viewer.setUpViewInWindow(100, 100, 800, 600);
    return viewer.run();
}
