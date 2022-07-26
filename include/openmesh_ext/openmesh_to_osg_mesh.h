#ifndef OPENMESH_TO_OSG_MESH_H_
#define OPENMESH_TO_OSG_MESH_H_

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include <osg/Geometry>
#include <osg/Array>

#include "openmesh_typedef.h"

void TriMeshToGeometry(TriMesh& triMesh, osg::ref_ptr<osg::Geometry>& osgTriMesh)
{
    triMesh.request_vertex_normals();
    triMesh.request_face_normals();
    triMesh.update_normals();

    osg::ref_ptr<osg::Vec3Array> vertices = new osg::Vec3Array(triMesh.n_vertices());
    osg::ref_ptr<osg::Vec3Array> normals = new osg::Vec3Array(triMesh.n_vertices());
    for (int i = 0; i < triMesh.n_vertices(); ++i)
    {
        auto& pt = triMesh.point(triMesh.vertex_handle(i));
        (*vertices)[i].set(pt[0], pt[1], pt[2]);
        auto& normal = triMesh.normal(triMesh.vertex_handle(i));
        (*normals)[i].set(normal[0], normal[1], normal[2]);
    }

    osg::ref_ptr<osg::DrawElementsUInt> indices = new osg::DrawElementsUInt(GL_TRIANGLES, triMesh.n_faces() * 3);
    for (int i = 0; i < triMesh.n_faces(); ++i)
    {
        int j = 0;
        auto& fh = triMesh.face_handle(i);
        for (auto iter = triMesh.fv_begin(fh); iter != triMesh.fv_end(fh); ++iter)
        {
            (*indices)[i * 3 + j] = iter->idx();
            j++;
        }
    }

    osgTriMesh->setDataVariance(osg::Object::DYNAMIC);
    osgTriMesh->setUseDisplayList(false);
    osgTriMesh->setUseVertexBufferObjects(true);
    osgTriMesh->setUseVertexArrayObject(true);

    osgTriMesh->setVertexArray(vertices.get());
    osgTriMesh->setNormalArray(normals.get());
    osgTriMesh->setNormalBinding(osg::Geometry::BIND_PER_VERTEX);
    osgTriMesh->addPrimitiveSet(indices.get());
}

bool loadTriangleMesh(const std::string& filename, osg::ref_ptr<osg::Geometry>& osgTriMesh)
{
    TriMesh triMesh;
    if (!OpenMesh::IO::read_mesh(triMesh, filename))
    {
        return false;
    }

    TriMeshToGeometry(triMesh, osgTriMesh);

    return true;
}


#endif
