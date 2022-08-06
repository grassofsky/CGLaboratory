#ifndef OPENMESH_UTILS_H_
#define OPENMESH_UTILS_H_

#include "openmesh_typedef.h"

void GetTriMeshBoundingBox(const TriMesh& trimesh, TriMesh::Point& low, TriMesh::Point& high)
{
    auto& pt_first = trimesh.point(trimesh.vertices().begin());
    low = pt_first;
    high = pt_first;
    for (auto it : trimesh.vertices())
    {
        const TriMesh::Point& point = trimesh.point(it);
        if (low[0] > point[0]) low[0] = point[0];
        if (low[1] > point[1]) low[1] = point[1];
        if (low[2] > point[2]) low[2] = point[2];

        if (high[0] < point[0]) high[0] = point[0];
        if (high[1] < point[1]) high[1] = point[1];
        if (high[2] < point[2]) high[2] = point[2];
    }
}

#endif