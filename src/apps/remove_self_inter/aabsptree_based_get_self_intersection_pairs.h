#include "openmesh_ext/openmesh_axis_aligned_bsptree.h"
#include "Mathematics/IntrSegment3Triangle3.h"
#include "Mathematics/IntrTriangle3Triangle3.h"

namespace cglab
{
namespace details
{
    void GetSharedVertices(std::array<TriMesh::VertexHandle, 3>& f0, std::array<TriMesh::VertexHandle, 3>& f1, std::vector<TriMesh::VertexHandle>& sharedVertices)
    {
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                if (f0[i] == f1[j])
                {
                    sharedVertices.push_back(f0[i]);
                }
            }
        }
    }
} // namespace details

std::vector<std::pair<int, int>> GetSelfIntersectionPairsAABspTreeBased(TriMesh& mesh)
{
    auto pt = mesh.point(mesh.vertices().begin());
    TriMesh::Point minPt = pt;
    TriMesh::Point maxPt = pt;
    for (auto vertex : mesh.vertices())
    {
        auto& ptTmp = mesh.point(vertex);
        if (minPt[0] > ptTmp[0]) minPt[0] = ptTmp[0];
        if (minPt[1] > ptTmp[1]) minPt[1] = ptTmp[1];
        if (minPt[2] > ptTmp[2]) minPt[2] = ptTmp[2];

        if (maxPt[0] < ptTmp[0]) maxPt[0] = ptTmp[0];
        if (maxPt[1] < ptTmp[1]) maxPt[1] = ptTmp[1];
        if (maxPt[2] < ptTmp[2]) maxPt[2] = ptTmp[2];
    }
    openmesh_ext::AABspTree tree;
    tree.GenerateBspTree(mesh, minPt, maxPt, 16, 16);

    const std::vector<openmesh_ext::AABspTreeNode *>& leafCollection = tree.GetLeafCollection();

    std::array<TriMesh::VertexHandle, 3> f0_v, f1_v;
    gte::TIQuery<float, gte::Triangle3<float>, gte::Triangle3<float>> queryTri3Tri3;
    gte::TIQuery<float, gte::Segment3<float>, gte::Triangle3<float>> querySegment3Tri3;
    std::set<std::pair<int, int>> interResult;
    for (int iNode = 0; iNode < leafCollection.size(); ++iNode)
    {
        auto& faceIdxs = leafCollection[iNode]->faces;
        if (faceIdxs.empty())
        {
            continue;
        }
        for (int iFace = 0; iFace < faceIdxs.size() - 1; ++iFace)
        {
            for (int jFace = iFace + 1; jFace < faceIdxs.size(); ++jFace)
            {
                auto f0v_iter = mesh.fv_iter(mesh.face_handle(faceIdxs[iFace]));
                f0_v[0] = *(f0v_iter++);
                f0_v[1] = *(f0v_iter++);
                f0_v[2] = *(f0v_iter++);

                auto f1v_iter = mesh.fv_iter(mesh.face_handle(faceIdxs[jFace]));
                f1_v[0] = *(f1v_iter++);
                f1_v[1] = *(f1v_iter++);
                f1_v[2] = *(f1v_iter++);

                std::vector<TriMesh::VertexHandle> sharedVertices;
                details::GetSharedVertices(f0_v, f1_v, sharedVertices);

                size_t sharedCount = sharedVertices.size();

                bool intersect = false;

                if (sharedCount == 3)
                {
                    intersect = true;
                }
                else if (sharedCount != 2)
                {
                    auto& f0_pt0 = mesh.point(f0_v[0]);
                    auto& f0_pt1 = mesh.point(f0_v[1]);
                    auto& f0_pt2 = mesh.point(f0_v[2]);

                    auto& f1_pt0 = mesh.point(f1_v[0]);
                    auto& f1_pt1 = mesh.point(f1_v[1]);
                    auto& f1_pt2 = mesh.point(f1_v[2]);

                    gte::Triangle3<float> tri0, tri1;
                    tri0.v[0] = { f0_pt0[0], f0_pt0[1], f0_pt0[2] };
                    tri0.v[1] = { f0_pt1[0], f0_pt1[1], f0_pt1[2] };
                    tri0.v[2] = { f0_pt2[0], f0_pt2[1], f0_pt2[2] };

                    tri1.v[0] = { f1_pt0[0], f1_pt0[1], f1_pt0[2] };
                    tri1.v[1] = { f1_pt1[0], f1_pt1[1], f1_pt1[2] };
                    tri1.v[2] = { f1_pt2[0], f1_pt2[1], f1_pt2[2] };

                    if (sharedCount == 0)
                    {
                        auto result = queryTri3Tri3(tri0, tri1);
                        intersect = result.intersect;
                    }
                    else if (sharedCount == 1)
                    {
                        std::vector<TriMesh::Point> f0_edge, f1_edge;
                        if (f0_v[0] != sharedVertices[0]) f0_edge.push_back(mesh.point(f0_v[0]));
                        if (f0_v[1] != sharedVertices[0]) f0_edge.push_back(mesh.point(f0_v[1]));
                        if (f0_v[2] != sharedVertices[0]) f0_edge.push_back(mesh.point(f0_v[2]));

                        if (f1_v[0] != sharedVertices[0]) f1_edge.push_back(mesh.point(f1_v[0]));
                        if (f1_v[1] != sharedVertices[0]) f1_edge.push_back(mesh.point(f1_v[1]));
                        if (f1_v[2] != sharedVertices[0]) f1_edge.push_back(mesh.point(f1_v[2]));

                        if (f0_edge.size() == 2 && f1_edge.size() == 2)
                        {
                            gte::Segment3<float> f0_seg, f1_seg;
                            f0_seg.p[0] = { f0_edge[0][0], f0_edge[0][1], f0_edge[0][2] };
                            f0_seg.p[1] = { f0_edge[1][0], f0_edge[1][1], f0_edge[1][2] };

                            f1_seg.p[0] = { f1_edge[0][0], f1_edge[0][1], f1_edge[0][2] };
                            f1_seg.p[1] = { f1_edge[1][0], f1_edge[1][1], f1_edge[1][2] };

                            auto result = querySegment3Tri3(f0_seg, tri1);
                            intersect = result.intersect;
                            result = querySegment3Tri3(f1_seg, tri0);
                            intersect = intersect || result.intersect;
                        }
                    }
                }

                if (intersect)
                {
                    std::pair<int, int> interPair;
                    if (faceIdxs[iFace] < faceIdxs[jFace])
                    {
                        interPair.first = faceIdxs[iFace];
                        interPair.second = faceIdxs[jFace];
                    }
                    else
                    {
                        interPair.first = faceIdxs[jFace];
                        interPair.second = faceIdxs[iFace];
                    }

                    interResult.insert(interPair);
                }
            }
        }
    }
    tree.DeleteBspTree();

    return std::vector<std::pair<int, int>>(interResult.begin(), interResult.end());
}
} // namespace cglab