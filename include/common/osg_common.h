#include <osg/Camera>
#include <osgText/Text>
#include <osgText/Font>

#include "colormap_define.h"

#ifndef OSG_COMMON_H_
#define OSG_COMMON_H_

osg::ref_ptr<osgText::Font> g_font = osgText::readFontFile("fonts/arial.ttf");

osg::Camera* CreateHUDCamera( double left, double right, double bottom, double top )
{
    osg::ref_ptr<osg::Camera> camera = new osg::Camera;
    camera->setReferenceFrame( osg::Transform::ABSOLUTE_RF );
    camera->setClearMask( GL_DEPTH_BUFFER_BIT );
    camera->setRenderOrder( osg::Camera::POST_RENDER );
    camera->setAllowEventFocus( false );
    camera->setProjectionMatrix( osg::Matrix::ortho2D(left, right, bottom, top) );
    camera->getOrCreateStateSet()->setMode( GL_LIGHTING, osg::StateAttribute::OFF );
    return camera.release();
}

osgText::Text* CreateText( const osg::Vec3& pos, const std::string& content, float size )
{
    osg::ref_ptr<osgText::Text> text = new osgText::Text;
    text->setDataVariance( osg::Object::DYNAMIC );
    text->setFont( g_font.get() );
    text->setCharacterSize( size );
    text->setAxisAlignment( osgText::TextBase::XY_PLANE );
    text->setPosition( pos );
    text->setText( content );
    return text.release();
}

// The color map's color changes from bottom to top
osg::Geometry* CreateColorMap(double left, double right, double bottom, double top)
{
    osg::ref_ptr<osg::Vec3Array> vertices = new osg::Vec3Array;

    double step = 4.0;

    for (double istart = bottom; istart < top; istart += step)
    {
        vertices->push_back(osg::Vec3(left, istart, 0));
        vertices->push_back(osg::Vec3(right, istart, 0));
    }

    double red, green, blue;
    osg::ref_ptr<osg::Vec3Array> colors = new osg::Vec3Array(vertices->size());
    for (int i = 0; i < colors->size(); ++i)
    {
        auto& pos = (*vertices)[i];
        InterplateColor((pos[1] - bottom)/(top - bottom), red, green, blue);
        (*colors)[i].set(red, green, blue);
    }

    osg::ref_ptr<osg::DrawElementsUInt> indices = new osg::DrawElementsUInt(GL_TRIANGLES, (vertices->size()  - 2) * 3);
    for (int i = 0; i < vertices->size() - 2; ++i)
    {
        (*indices)[i * 3] = i; 
        (*indices)[i * 3 + 1] = i + 1;
        (*indices)[i * 3 + 2] = i + 2;
    }

    osg::ref_ptr<osg::Vec3Array> normal = new osg::Vec3Array();
    normal->push_back(osg::Vec3(0.0, 0.0, 1.0));

    osg::ref_ptr<osg::Geometry> geom = new osg::Geometry;
    geom->setVertexArray(vertices.get());
    geom->setColorArray(colors.get());
    geom->setNormalArray(normal);
    geom->setNormalBinding(osg::Geometry::BIND_OVERALL);
    geom->setColorBinding(osg::Geometry::BIND_PER_VERTEX);
    geom->addPrimitiveSet(indices.get());

    return geom.release();
}

#endif