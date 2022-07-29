/// see: https://learnopengl-cn.github.io/01%20Getting%20started/06%20Textures/

#include "config.h"

#include <sstream>

#include <osg/Program>
#include <osg/Geometry>
#include <osg/Geode>
#include <osg/Group>
#include <osg/Image>
#include <osg/Texture2D>
#include <osgViewer/Viewer>


#define STB_IMAGE_IMPLEMENTATION
#include "stb_image/stb_image.h"

// in class drawable
// enum AttributeTypes
//{
//    VERTICES = 0,
//    WEIGHTS = 1,
//    NORMALS = 2,
//    COLORS = 3,
//    SECONDARY_COLORS = 4,
//    FOG_COORDS = 5,
//    ATTRIBUTE_6 = 6,
//    ATTRIBUTE_7 = 7,
//    TEXTURE_COORDS = 8,
//    TEXTURE_COORDS_0 = TEXTURE_COORDS,
//    TEXTURE_COORDS_1 = TEXTURE_COORDS_0 + 1,
//    TEXTURE_COORDS_2 = TEXTURE_COORDS_0 + 2,
//    TEXTURE_COORDS_3 = TEXTURE_COORDS_0 + 3,
//    TEXTURE_COORDS_4 = TEXTURE_COORDS_0 + 4,
//    TEXTURE_COORDS_5 = TEXTURE_COORDS_0 + 5,
//    TEXTURE_COORDS_6 = TEXTURE_COORDS_0 + 6,
//    TEXTURE_COORDS_7 = TEXTURE_COORDS_0 + 7
//    // only eight texture coord examples provided here, but underlying code can handle any no of texture units,
//    // simply co them as (TEXTURE_COORDS_0+unit).
//};
const char* vertexShaderSource = "#version 330 core\n"
"layout (location = 0) in vec4 aPos;\n"
"layout (location = 3) in vec4 aColor;\n"
"layout (location = 8) in vec2 aTexCoord;\n"
"out vec4 v2fColor;\n"
"out vec2 v2fTexCoord;\n"
"void main()\n"
"{\n"
"   gl_Position = aPos;\n"
"   v2fColor = aColor;\n"
"   v2fTexCoord = aTexCoord;\n"
"}\0";
const char* fragmentShaderSource = "#version 330 core\n"
"out vec4 FragColor;\n"
"in vec4 v2fColor;\n"
"in vec2 v2fTexCoord;\n"
"uniform vec4 ourColorFactor;\n"
"uniform sampler2D texture1;\n"
"void main()\n"
"{\n"
"   FragColor = v2fColor * ourColorFactor * texture(texture1, v2fTexCoord);\n"
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

    osg::ref_ptr<osg::Vec2Array> texcoords = new osg::Vec2Array(3);
    (*texcoords)[0].set(0.0, 0.0);
    (*texcoords)[1].set(1.0, 0.0);
    (*texcoords)[2].set(0.5, 1.0);

    osg::ref_ptr<osg::Geometry> geoTriangle = new osg::Geometry;
    geoTriangle->setDataVariance(osg::Object::DYNAMIC);
    geoTriangle->setUseDisplayList(false);
    geoTriangle->setUseVertexBufferObjects(true);
    geoTriangle->setUseVertexArrayObject(true);
    geoTriangle->setVertexArray(vertices.get());
    geoTriangle->addPrimitiveSet(indices.get());
    geoTriangle->setColorArray(colors.get());
    geoTriangle->setTexCoordArray(0, texcoords.get());

    geoTriangle->setColorBinding(osg::Geometry::BIND_PER_VERTEX);

    osg::ref_ptr<osg::Geode> geode = new osg::Geode;
    geode->addDrawable(geoTriangle);

    osg::ref_ptr<osg::Program> program = new osg::Program();
    program->addShader(new osg::Shader(osg::Shader::VERTEX, vertexShaderSource));
    program->addShader(new osg::Shader(osg::Shader::FRAGMENT, fragmentShaderSource));

    // set program
    geode->getOrCreateStateSet()->setAttributeAndModes(program, osg::StateAttribute::ON);
    // Set uniform value
    geode->getOrCreateStateSet()->addUniform(new osg::Uniform("ourColorFactor", osg::Vec4f(1.0, 1.0, 1.0, 1.0)));
    

    // load texture image
    int width, height, nrChannels;
    std::stringstream ss;
    ss << TEXTURE_PATH << "container.jpg";
    unsigned char* pImageData = stbi_load(ss.str().c_str(), &width, &height, &nrChannels, 0);
    osg::ref_ptr<osg::Image> textureImage = new osg::Image;
    textureImage->setImage(width, height, 1, nrChannels == 3 ? GL_RGB : GL_RGBA, nrChannels == 3 ? GL_RGB : GL_RGBA, GL_UNSIGNED_BYTE, pImageData, osg::Image::NO_DELETE);
    osg::ref_ptr<osg::Texture2D> texture = new osg::Texture2D;
    texture->setImage(textureImage.get());

    geode->getOrCreateStateSet()->setTextureAttributeAndModes(0, texture.get(), osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);

    osg::ref_ptr<osg::Group> root = new osg::Group;
    root->addChild(geode);

    // transform by camera
    osgViewer::Viewer viewer;
    viewer.setSceneData(root);
    viewer.setUpViewInWindow(100, 100, 800, 600);
    return viewer.run();
}
