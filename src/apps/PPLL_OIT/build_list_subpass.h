#ifndef BUILD_LIST_SUBPASS_H_
#define BUILD_LIST_SUBPASS_H_

#include <osg/Group>
#include <osg/Texture2D>
#include <osg/BindImageTexture>

class BuildListSubpass : public osg::Group
{
public:
    BuildListSubpass();
    ~BuildListSubpass();

private:
    CreateScene_();

    unsigned int view_port_width_;
    unsigned int view_port_height_;
    osg::ref_ptr<osg::Texture2D> head_pointer_texture_;
    osg::ref_ptr<osg::BindImageTexture> head_pointer_image_bind_;
};

#endif
