

#include <osg/Geometry>
#include <osg/Program>
#include <osg/ShapeDrawable>
#include <osg/Geode>
#include <osg/LineWidth>
#include <osg/Camera>
#include <osg/Texture2D>
#include <osg/PolygonMode>
#include <osgDB/ReadFile>
#include <osgUtil/TangentSpaceGenerator>
#include <osgViewer/Viewer>

std::string strDataPath = "E:/workplace/github/OpenSceneGraph-Data/";

osg::Geode* createScreenQuad(float width, float height,
    float scale)
{
    osg::Geometry* geom = osg::createTexturedQuadGeometry(
        osg::Vec3(), osg::Vec3(width, 0.0f, 0.0f),
        osg::Vec3(0.0f, height, 0.0f),
        0.0f, 0.0f, width * scale, height * scale);
    geom->setUseVertexArrayObject(true);
    geom->setUseDisplayList(false);
    osg::ref_ptr<osg::Geode> quad = new osg::Geode;
    quad->addDrawable(geom);
    int values = osg::StateAttribute::OFF |
        osg::StateAttribute::PROTECTED;
    quad->getOrCreateStateSet()->setAttribute(
        new osg::PolygonMode(osg::PolygonMode::FRONT_AND_BACK,
            osg::PolygonMode::FILL), values);
    quad->getOrCreateStateSet()->setMode(GL_LIGHTING, values);
    return quad.release();
}

osg::Camera* createRTTCamera(osg::Camera::BufferComponent
    buffer, osg::Texture* tex, bool isAbsolute)
{
    osg::ref_ptr<osg::Camera> camera = new osg::Camera;
    camera->setClearColor(osg::Vec4());
    camera->setClearMask(
        GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    camera->setRenderTargetImplementation(
        osg::Camera::FRAME_BUFFER_OBJECT);
    camera->setRenderOrder(osg::Camera::PRE_RENDER);
    if (tex)
    {
        tex->setFilter(osg::Texture2D::MIN_FILTER,
            osg::Texture2D::LINEAR);
        tex->setFilter(osg::Texture2D::MAG_FILTER,
            osg::Texture2D::LINEAR);
        camera->setViewport(0, 0, tex->getTextureWidth(),
            tex->getTextureHeight());
        camera->attach(buffer, tex);
    }
    if (isAbsolute)
    {
        camera->setReferenceFrame(osg::Transform::ABSOLUTE_RF);
        camera->setProjectionMatrix(osg::Matrix::ortho2D(
            0.0, 1.0, 0.0, 1.0));
        camera->setViewMatrix(osg::Matrix::identity());
        camera->addChild(createScreenQuad(1.0f, 1.0f, 1.0f));
    }
    return camera.release();
}

static const char* vertSource = {
    "#version 330 \n"
    "layout(location = 0) in vec4 "
    "layout(location = 6) in vec3 tangent;\n"
    "layout(location = 7) in vec3 binormal;\n"
    "out vec3 lightDir;\n"
    "void main()\n"
    "{\n"
    "    vec3 normal = normalize(gl_NormalMatrix * gl_Normal);\n"
    "    mat3 rotation = mat3(tangent, binormal, normal);\n"
    "    vec4 vertexInEye = gl_ModelViewMatrix * gl_Vertex;\n"
    "    lightDir = vec3(gl_LightSource[0].position.xyz - vertexInEye.xyz);\n"
    "    lightDir = normalize(rotation * normalize(lightDir));\n"
    "    gl_Position = ftransform();\n"
    "    gl_TexCoord[0] = gl_MultiTexCoord0;\n"
    "}\n"
};

static const char* fragSource = {
    "#version 330\n"
    "uniform sampler2D colorTex;\n"
    "uniform sampler2D normalTex;\n"
    "in vec3 lightDir;\n"
    "void main (void)\n"
    "{\n"
    "    vec4 base = texture2D(colorTex, gl_TexCoord[0].xy);\n"
    "    vec3 bump = texture2D(normalTex, gl_TexCoord[0].xy).xyz;\n"
    "    bump = normalize(bump * 2.0 - 1.0);\n"

    "    float lambert = max(dot(bump, lightDir), 0.0);\n"
    "    if (lambert > 0.0)\n"
    "    {\n"
    "        gl_FragColor = base * gl_LightSource[0].diffuse * lambert;\n"
    "        gl_FragColor += gl_LightSource[0].specular * pow(lambert, 2.0);\n"
    "    }\n"
    "    gl_FragColor += gl_LightSource[0].ambient;\n"
    "}\n"
};

void generateTangentArray(osg::Geometry* geom)
{
    osg::ref_ptr<osgUtil::TangentSpaceGenerator> tsg =
        new osgUtil::TangentSpaceGenerator;
    tsg->generate(geom);
    geom->setVertexAttribArray(6, tsg->getTangentArray());
    geom->setVertexAttribBinding(
        6, osg::Geometry::BIND_PER_VERTEX);
    geom->setVertexAttribArray(7, tsg->getBinormalArray());
    geom->setVertexAttribBinding(
        7, osg::Geometry::BIND_PER_VERTEX);
}

class ComputeTangentVisitor : public osg::NodeVisitor
{
public:
    void apply(osg::Node& node) { traverse(node); }

    void apply(osg::Geode& node)
    {
        for (unsigned int i = 0; i < node.getNumDrawables(); ++i)
        {
            osg::Geometry* geom = dynamic_cast<osg::Geometry*>(node.getDrawable(i));
            if (geom) generateTangentArray(geom);
        }

        traverse(node);
    }
};

int main()
{
    osg::ref_ptr<osg::Node> scene = osgDB::readNodeFile(strDataPath + "skydome.osgt");
    ComputeTangentVisitor ctv;
    ctv.setTraversalMode(osg::NodeVisitor::TRAVERSE_ALL_CHILDREN);
    scene->accept(ctv);

    osg::ref_ptr<osg::Program> program = new osg::Program;
    program->addShader(new osg::Shader(osg::Shader::VERTEX,
        vertSource));
    program->addShader(new osg::Shader(osg::Shader::FRAGMENT,
        fragSource));
    program->addBindAttribLocation("tangent", 6);
    program->addBindAttribLocation("binormal", 7);

    osg::ref_ptr<osg::Texture2D> colorTex = new osg::Texture2D;
    colorTex->setImage(osgDB::readImageFile(
        strDataPath + "Images/whitemetal_diffuse.jpg"));
    osg::ref_ptr<osg::Texture2D> normalTex =
        new osg::Texture2D;
    normalTex->setImage(osgDB::readImageFile(
        strDataPath + "Images/whitemetal_normal.jpg"));
    osg::StateSet* stateset = scene->getOrCreateStateSet();
    stateset->addUniform(new osg::Uniform("colorTex", 0));
    stateset->addUniform(new osg::Uniform("normalTex", 1));
    stateset->setAttributeAndModes(program.get());
    osg::StateAttribute::GLModeValue value =
        osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE;
    stateset->setTextureAttributeAndModes(0, colorTex.get(),
        value);
    stateset->setTextureAttributeAndModes(1, normalTex.get(),
        value);

    osgViewer::Viewer viewer;
    viewer.setSceneData(scene.get());
    return viewer.run();
}