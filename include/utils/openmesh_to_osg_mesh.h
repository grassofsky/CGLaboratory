#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include <osg/Shape>
#include <osg/Array>

typedef OpenMesh::TriMesh_ArrayKernelT<> TriMesh;

bool loadTriangleMesh(const std::string& filename, osg::ref_ptr<osg::TriangleMesh>& osgTriMesh)
{
    TriMesh triMesh;
    if (!OpenMesh::IO::read_mesh(triMesh, filename))
    {
        return false;
    }

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

    osg::ref_ptr<osg::UIntArray> indices = new osg::UIntArray(triMesh.n_faces() * 3);
    for (int i = 0; i < triMesh.n_faces(); ++i)
    {
        int j = 0;
        for (auto iter = triMesh.fv_begin(triMesh.face_handle(i)); iter != triMesh.fv_end(); ++iter)
        {
            (*indices)[i * 3 + j] = iter->idx();
            j++;
        }
    }

    osgTriMesh->setVertices(vertices.get());
    osgTriMesh->setIndices(indices.get());
}

