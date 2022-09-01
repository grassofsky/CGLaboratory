#include "build_list_subpass.h"

BuildListSubpass::BuildListSubpass()
: view_port_height_(100)
, view_port_width_(100)
{
}

BuildListSubpass::~BuildListSubpass()
{
}

void BuildListSubpass::CreateScene_()
{
    osg::ref_ptr<osg::Image> image = new osg::Image;
    image->setImage(view_port_width_, view_port_height_, 1, GL_R32UI, GL_RED_INTEGER, GL_UNSIGNED_INT, nullptr, osg::Image::USE_MALLOC_FREE)
    head_pointer_texture_ = new osg::Texture2D();

}

