#include <iostream>
#include <fstream>
#include <unordered_set>
#include <unordered_map>

#include "openmesh_ext/openmesh_utils.h"
#include "Mathematics/DistPointTriangle.h"

#include "config.h"

/// Fast computing adaptively sampled distance field on GPU - 2011

typedef std::pair<int, TriMesh::FaceHandle> TypeVoxelFacePair;
struct TypeVoxelFacePairHashFunc
{
    size_t operator()(const TypeVoxelFacePair& pair) const
    {
        return std::hash<int>()(pair.first) ^ std::hash<int>()(pair.second.idx());
    }
};

const unsigned int knmax = 16;

double GetPointToFaceDistance(TriMesh::Point& pt, TriMesh& trimesh2, TriMesh::FaceHandle fh)
{
    gte::DCPPoint3Triangle3<double> distanceQuery;
    gte::Vector<3, double> point = { pt[0], pt[1], pt[2] };

    TriMesh::Point triPt[3];
    auto fviter = trimesh2.fv_iter(fh);
    triPt[0] = trimesh2.point(fviter++);
    triPt[1] = trimesh2.point(fviter++);
    triPt[2] = trimesh2.point(fviter++);

    gte::Triangle<3, double> triangle;
    triangle.v[0] = { triPt[0][0], triPt[0][1], triPt[0][2] };
    triangle.v[1] = { triPt[1][0], triPt[1][1], triPt[1][2] };
    triangle.v[2] = { triPt[2][0], triPt[2][1], triPt[2][2] };

    return distanceQuery(point, triangle).distance;
}

TriMesh::Point StaticCastToInt(const TriMesh::Point& input)
{
    TriMesh::Point out;
    out[0] = static_cast<int>(input[0]);
    out[1] = static_cast<int>(input[1]);
    out[2] = static_cast<int>(input[2]);
    return out;
}

// 两个包围盒，low的最大值要小于high的最小值
bool CheckIfAABBIntersect(const TriMesh::Point& low1, const TriMesh::Point& high1, const TriMesh::Point& low2, const TriMesh::Point& high2)
{
    return std::max(low1[0], low2[0]) <= std::min(high1[0], high2[0]) &&
        std::max(low1[1], low2[1]) <= std::min(high1[1], high2[1]) &&
        std::max(low1[2], low2[2]) <= std::min(high1[2], high2[2]);
}

bool CheckIfAABBIntersect(const double& low1x, const double& low1y, const double& low1z, 
    const double& high1x, const double& high1y, const double& high1z,
    const double& low2x, const double& low2y, const double& low2z,
    const double& high2x, const double& high2y, const double& high2z)
{
    return std::max(low1x, low2x) <= std::min(high1x, high2x) &&
        std::max(low1y, low2y) <= std::min(high1y, high2y) &&
        std::max(low1z, low2z) <= std::min(high1z, high2z);
}

bool CheckIfVoxelTriangleIntersect(const TriMesh::Point& low, const TriMesh::Point& high, const TriMesh& trimesh, TriMesh::FaceHandle fh)
{
    auto fviter = trimesh.cfv_iter(fh);
    auto& v0 = trimesh.point(fviter++);
    auto& v1 = trimesh.point(fviter++);
    auto& v2 = trimesh.point(fviter++);

    auto delta = high - low;
    TriMesh::Point c;
    auto& normal = trimesh.normal(fh);
    c[0] = normal[0] > 0 ? delta[0] : 0;
    c[1] = normal[1] > 0 ? delta[1] : 0;
    c[2] = normal[2] > 0 ? delta[2] : 0;

    // check if voxel intersect with plane of triangle
    auto d1 = normal.dot(c - v0);
    auto d2 = normal.dot((delta - c) - v0);
    auto ndotp = normal.dot(low);
    if (!((ndotp + d1) * (ndotp + d2) <= 0))
    {
        return false;
    }

    // check projection intersection
    auto e0 = v1 - v0;
    auto e1 = v2 - v1;
    auto e2 = v0 - v2;

    // xy projection
    auto sign_xy = (normal[2] >= 0 ? 1 : -1);
    auto n_xy_e0 = TriMesh::Point(-e0[1], e0[0], 0) * sign_xy;
    auto n_xy_e1 = TriMesh::Point(-e1[1], e1[0], 0) * sign_xy;
    auto n_xy_e2 = TriMesh::Point(-e2[1], e2[0], 0) * sign_xy;

    auto d_xy_e0 = -n_xy_e0.dot(TriMesh::Point(v0[0], v0[1], 0)) + std::max(0.f, delta[0] * n_xy_e0[0]) + std::max(0.f, delta[1] * n_xy_e0[1]);
    auto d_xy_e1 = -n_xy_e1.dot(TriMesh::Point(v1[0], v1[1], 0)) + std::max(0.f, delta[0] * n_xy_e1[0]) + std::max(0.f, delta[1] * n_xy_e1[1]);
    auto d_xy_e2 = -n_xy_e2.dot(TriMesh::Point(v2[0], v2[1], 0)) + std::max(0.f, delta[0] * n_xy_e2[0]) + std::max(0.f, delta[1] * n_xy_e2[1]);

    bool xy_ok = (n_xy_e0.dot(low) + d_xy_e0 >= 0) && (n_xy_e1.dot(low) + d_xy_e1 >= 0) && (n_xy_e2.dot(low) + d_xy_e2 >= 0);

    if (xy_ok)
        return true;

    // yz projection
    auto sign_yz = (normal[0] >= 0 ? 1 : -1);
    auto n_yz_e0 = TriMesh::Point(-e0[2], e0[1], 0) * sign_yz;
    auto n_yz_e1 = TriMesh::Point(-e1[2], e1[1], 0) * sign_yz;
    auto n_yz_e2 = TriMesh::Point(-e2[2], e2[1], 0) * sign_yz;

    auto d_yz_e0 = -n_yz_e0.dot(TriMesh::Point(v0[1], v0[2], 0)) + std::max(0.f, delta[1] * n_yz_e0[1]) + std::max(0.f, delta[2] * n_yz_e0[2]);
    auto d_yz_e1 = -n_yz_e1.dot(TriMesh::Point(v1[1], v1[2], 0)) + std::max(0.f, delta[1] * n_yz_e1[1]) + std::max(0.f, delta[2] * n_yz_e1[2]);
    auto d_yz_e2 = -n_yz_e2.dot(TriMesh::Point(v2[1], v2[2], 0)) + std::max(0.f, delta[1] * n_yz_e2[1]) + std::max(0.f, delta[2] * n_yz_e2[2]);

    bool yz_ok = (n_yz_e0.dot(low) + d_yz_e0 >= 0) && (n_yz_e1.dot(low) + d_yz_e1 >= 0) && (n_yz_e2.dot(low) + d_yz_e2 >= 0);
    if (yz_ok)
        return true;


    // zx projection
    auto sign_zx = (normal[1] >= 0 ? 1 : -1);
    auto n_zx_e0 = TriMesh::Point(-e0[0], e0[2], 0) * sign_zx;
    auto n_zx_e1 = TriMesh::Point(-e1[0], e1[2], 0) * sign_zx;
    auto n_zx_e2 = TriMesh::Point(-e2[0], e2[2], 0) * sign_zx;

    auto d_zx_e0 = -n_zx_e0.dot(TriMesh::Point(v0[2], v0[0], 0)) + std::max(0.f, delta[0] * n_zx_e0[0]) + std::max(0.f, delta[2] * n_zx_e0[2]);
    auto d_zx_e1 = -n_zx_e1.dot(TriMesh::Point(v1[2], v1[0], 0)) + std::max(0.f, delta[0] * n_zx_e1[0]) + std::max(0.f, delta[2] * n_zx_e1[2]);
    auto d_zx_e2 = -n_zx_e2.dot(TriMesh::Point(v2[2], v2[0], 0)) + std::max(0.f, delta[0] * n_zx_e2[0]) + std::max(0.f, delta[2] * n_zx_e2[2]);

    bool zx_ok = (n_zx_e0.dot(low) + d_zx_e0 >= 0) && (n_zx_e1.dot(low) + d_zx_e1 >= 0) && (n_zx_e2.dot(low) + d_zx_e2 >= 0);

    //return xy_ok || yz_ok || zx_ok;
    return zx_ok;
}


int main()
{
    // load trimesh and calculate distance
    TriMesh trimesh;
    OpenMesh::IO::read_mesh(trimesh, std::string(MESHES_PATH) + "rabit.off");
    trimesh.request_vertex_normals();
    trimesh.request_face_normals();
    trimesh.update_normals();

    auto start = std::chrono::steady_clock::now();

    // Init aabb infos
    OpenMesh::FPropHandleT<std::pair<TriMesh::Point, TriMesh::Point>> faceAABBProperty;
    trimesh.add_property(faceAABBProperty);
    for (auto fh : trimesh.faces())
    {
        auto fviter = trimesh.fv_iter(fh);
        auto& pt0 = trimesh.point(fviter++);
        auto& pt1 = trimesh.point(fviter++);
        auto& pt2 = trimesh.point(fviter++);
        std::pair<TriMesh::Point, TriMesh::Point> aabb;
        aabb.first = pt0;
        aabb.second = pt0;

        if (aabb.first[0] > pt1[0]) aabb.first[0] = pt1[0];
        if (aabb.first[1] > pt1[1]) aabb.first[1] = pt1[1];
        if (aabb.first[2] > pt1[2]) aabb.first[2] = pt1[2];
        if (aabb.first[0] > pt2[0]) aabb.first[0] = pt2[0];
        if (aabb.first[1] > pt2[1]) aabb.first[1] = pt2[1];
        if (aabb.first[2] > pt2[2]) aabb.first[2] = pt2[2];

        if (aabb.second[0] < pt1[0]) aabb.second[0] = pt1[0];
        if (aabb.second[1] < pt1[1]) aabb.second[1] = pt1[1];
        if (aabb.second[2] < pt1[2]) aabb.second[2] = pt1[2];
        if (aabb.second[0] < pt2[0]) aabb.second[0] = pt2[0];
        if (aabb.second[1] < pt2[1]) aabb.second[1] = pt2[1];
        if (aabb.second[2] < pt2[2]) aabb.second[2] = pt2[2];


        trimesh.property(faceAABBProperty, fh) = aabb;
    }

    // The corner point of volume
    TriMesh::Point low, high;
    GetTriMeshBoundingBox(trimesh, low, high);

    TriMesh::Point offset(0.1, 0.1, 0.1);
    low = low - offset;
    high = high + offset;

    int dim[] = { 512,512,512 };
    TriMesh::Point spacing = (high - low) / (TriMesh::Point(dim[0], dim[1], dim[2]));

    TriMesh::Point extLow = low - spacing * 0.5;

    int slicesize = dim[0] * dim[1];
    int volumesize = slicesize * dim[2];
    float* sdf = new float[volumesize];
    memset(sdf, 0, volumesize * sizeof(float));

    // loop all triangles


    std::unordered_map<int, std::vector<TriMesh::FaceHandle>> setIntersections;
    int nfaces = trimesh.n_faces();
//#pragma omp parallel for num_threads(8)
    for (int i = 0; i < nfaces; ++i)
    {
        std::pair<TriMesh::Point, TriMesh::Point>& aabb = trimesh.property(faceAABBProperty, trimesh.face_handle(i));

        // 根据boundingbox获取和aabb相交的voxel
        auto voxelLow = StaticCastToInt((aabb.first - extLow)/spacing) * spacing + extLow;
        auto voxelHigh = (StaticCastToInt((aabb.second - extLow) / spacing) + TriMesh::Point(1, 1, 1)) * spacing + extLow;

        for (float voxelx = voxelLow[0]; voxelx < voxelHigh[0]; voxelx += spacing[0])
        {
            for (float voxely = voxelLow[1]; voxely < voxelHigh[1]; voxely += spacing[1])
            {
                for (float voxelz = voxelLow[2]; voxelz < voxelHigh[2]; voxelz += spacing[2])
                {
                    TriMesh::Point pt000(voxelx, voxely, voxelz);

                    if (CheckIfAABBIntersect(pt000, pt000 + spacing, aabb.first, aabb.second))
                    {
                        //如果voxel和三角形的aabb相交，进一步判断voxel是否和三角形相交
                        //if (CheckIfVoxelTriangleIntersect(pt000, pt000 + spacing, trimesh, trimesh.face_handle(i)))
                        {
                            int xidx = static_cast<int>((pt000[0] - low[0]) / spacing[0]);
                            if (xidx < 0) continue;
                            else if (xidx > dim[0] - 1) continue;
                            int yidx = static_cast<int>((pt000[1] - low[1]) / spacing[1]);
                            if (yidx < 0) continue;
                            else if (yidx > dim[1] - 1) continue;
                            int zidx = static_cast<int>((pt000[2] - low[2]) / spacing[2]);
                            if (zidx < 0) continue;
                            else if (zidx > dim[2] - 1) continue;
                            int volumeidx = zidx * slicesize + yidx * dim[0] + xidx;
                            setIntersections[volumeidx].push_back(trimesh.face_handle(i));
                            //std::wcout << "True";
                            // 如果voxel和三角形相交，计算voxel的每个顶点到三角形的距离
                            //TriMesh::Point voxelhigh = pt000 + spacing;
                            //TriMesh::Point pt100(voxelhigh[0], voxely, voxelz);
                            //TriMesh::Point pt010(voxelx, voxelhigh[1], voxelz);
                            //TriMesh::Point pt001(voxelx, voxely, voxelhigh[2]);
                            //TriMesh::Point pt110(voxelhigh[0], voxelhigh[1], voxelz);
                            //TriMesh::Point pt101(voxelhigh[0], voxely, voxelhigh[2]);
                            //TriMesh::Point pt011(voxelx, voxelhigh[1], voxelhigh[2]);
                            //TriMesh::Point pt111(voxelhigh[0], voxelhigh[1], voxelhigh[2]);

                            //TriMesh::Point* pts[] = { &pt000, &pt100, &pt010, &pt001, &pt110, &pt101, &pt011, &pt111 };
                            //for (int ipt = 0; ipt < 8; ++ipt)
                            //{
                            //    int xidx = static_cast<int>(((*pts[ipt])[0] - low[0]) / spacing[0]);
                            //    if (xidx < 0) continue;
                            //    else if (xidx > dim[0] - 1) continue;
                            //    int yidx = static_cast<int>(((*pts[ipt])[1] - low[1]) / spacing[1]);
                            //    if (yidx < 0) continue;
                            //    else if (yidx > dim[1] - 1) continue;
                            //    int zidx = static_cast<int>(((*pts[ipt])[2] - low[2]) / spacing[2]);
                            //    if (zidx < 0) continue;
                            //    else if (zidx > dim[2] - 1) continue;
                            //    int volumeidx = zidx * slicesize + yidx * dim[0] + xidx;
                            //    float value = GetPointToFaceDistance(*pts[ipt], trimesh, trimesh.face_handle(i));
                            //    if (sdf[volumeidx] == 0.0f)
                            //    {
                            //        sdf[volumeidx] = value;
                            //    }
                            //    else if (value < sdf[volumeidx])
                            //    {
                            //        sdf[volumeidx] = value;
                            //    }
                            //}
                        }
                    }
                }
            }
        }

    }

    auto end = std::chrono::steady_clock::now();

    std::cout << "Generate volume cost: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        << " ms" << std::endl;


    std::ofstream ofile("D:/rabit.raw", std::ios::binary);
    ofile.write((char*)sdf, volumesize * sizeof(float));
    ofile.close();

    delete[] sdf;

    return 0;
}
