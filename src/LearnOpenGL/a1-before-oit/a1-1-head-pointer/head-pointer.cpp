#include <osg/Program>
#include <osg/Geometry>
#include <osg/Geode>
#include <osg/Group>
#include <osg/BufferIndexBinding>
#include <osg/BindImageTexture>
#include <osg/Texture2D>
#include <osgViewer/Viewer>

#include "config.h"
#include "openmesh_ext/openmesh_to_osg_mesh.h"

#include "head-pointer_vert.h"
#include "head-pointer_frag.h"

class ResetAtomicCounter : public osg::StateAttributeCallback
{
public:
    virtual void operator () (osg::StateAttribute* sa, osg::NodeVisitor*)
    {
        osg::AtomicCounterBufferBinding* acbb = dynamic_cast<osg::AtomicCounterBufferBinding*>(sa);
        if (acbb)
        {
            osg::BufferData* acbd = acbb->getBufferData();
            if (acbd)
            {
                acbd->dirty();
            }
        }
    }
};

class RetriveHeadPointerTexture : public osg::StateAttributeCallback
{
public:
    virtual void operator() (osg::StateAttribute* sa, osg::NodeVisitor*)
    {
        osg::BindImageTexture* bit = dynamic_cast<osg::BindImageTexture*>(sa);
        if (bit)
        {
            osg::Texture2D* tex = dynamic_cast<osg::Texture2D*>(bit->getTexture());
            if (tex)
            {
                // TODO
            }
        }
    }
};

int main()
{
    osg::ref_ptr<osg::Geometry> geoMesh = new osg::Geometry();
    if (!loadTriangleMesh(std::string(MESHES_PATH) + "sphere.stl", geoMesh))
    {
        return -1;
    }


    osg::ref_ptr<osg::Geode> geode = new osg::Geode;
    geode->addDrawable(geoMesh.get());

    osg::ref_ptr<osg::Program> program = new osg::Program();
    program->addShader(new osg::Shader(osg::Shader::VERTEX, vertexShaderSource));
    program->addShader(new osg::Shader(osg::Shader::FRAGMENT, fragmentShaderSource));
    geode->getOrCreateStateSet()->setAttributeAndModes(program, osg::StateAttribute::ON);

#pragma region Atomic
    osg::ref_ptr<osg::UIntArray> atomicCounterArray = new osg::UIntArray;
    atomicCounterArray->push_back(0);
    osg::ref_ptr<osg::AtomicCounterBufferObject> atomicCounterBO = new osg::AtomicCounterBufferObject;
    atomicCounterBO->setUsage(GL_STREAM_COPY);
    atomicCounterArray->setBufferObject(atomicCounterBO);

    osg::ref_ptr<osg::AtomicCounterBufferBinding> atomicCounterBB = new osg::AtomicCounterBufferBinding(0, atomicCounterArray.get(), 0, sizeof(GLuint));
    atomicCounterBB->setUpdateCallback(new ResetAtomicCounter);
    geode->getOrCreateStateSet()->setAttributeAndModes(atomicCounterBB);
#pragma endregion Atomic

#pragma region HeadPointer
    osg::ref_ptr<osg::Texture2D> headPointerTexture = new osg::Texture2D();
    const unsigned int kheadPointerTextureSize = 2048;
    headPointerTexture->setTextureSize(kheadPointerTextureSize, kheadPointerTextureSize);
    headPointerTexture->setFilter(osg::Texture2D::MIN_FILTER, osg::Texture2D::LINEAR);
    headPointerTexture->setFilter(osg::Texture2D::MAG_FILTER, osg::Texture2D::LINEAR);
    headPointerTexture->setInternalFormat(GL_R32UI);
    headPointerTexture->setSourceFormat(GL_RED_INTEGER_EXT);
    headPointerTexture->setSourceType(GL_UNSIGNED_INT);
    osg::ref_ptr<osg::BindImageTexture> headPointerBindImage = new osg::BindImageTexture(1,
        headPointerTexture,
        osg::BindImageTexture::READ_WRITE,
        GL_R32UI,
        0,
        false,
        0);
    geode->getOrCreateStateSet()->setAttributeAndModes(headPointerBindImage);
#pragma endregion HeadPointer

    osg::ref_ptr<osg::Group> root = new osg::Group;
    root->addChild(geode);

    osgViewer::Viewer viewer;
    viewer.setSceneData(root);
    viewer.setUpViewInWindow(100, 100, 800, 600);
    viewer.getCamera()->getGraphicsContext()->getState()->setUseModelViewAndProjectionUniforms(true);

    return viewer.run();
}
