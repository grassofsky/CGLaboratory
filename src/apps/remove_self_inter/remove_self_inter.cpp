#include <iostream>
#include <fstream>
#include <chrono>

#include "box_inter_based_get_self_intersection_pairs.h"
#include "aabsptree_based_get_self_intersection_pairs.h"

#include "config.h"

const unsigned int knmax = 16;

void WriteOutInterPairs(const std::string& outfilename, std::vector<std::pair<int,int>>& vecSelfInterPair, TriMesh& trimesh)
{
    std::ofstream ofile(outfilename);
    ofile << "solid\n";
    for (int i = 0; i < vecSelfInterPair.size(); ++i)
    {
        auto face1 = trimesh.face_handle(vecSelfInterPair[i].first);
        auto face2 = trimesh.face_handle(vecSelfInterPair[i].second);

        TriMesh::FaceHandle faces[] = { face1, face2 };
        for (int j = 0; j < 2; ++j)
        {
            auto normal = trimesh.normal(faces[j]);
            TriMesh::Point triPt[3];
            auto fviter = trimesh.fv_iter(faces[j]);
            triPt[0] = trimesh.point(fviter++);
            triPt[1] = trimesh.point(fviter++);
            triPt[2] = trimesh.point(fviter++);

            ofile << "facet normal " << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
            ofile << "outer loop\n";
            for (int k = 0; k < 3; ++k)
                ofile << "vertex " << triPt[k][0] << " " << triPt[k][1] << " " << triPt[k][2] << "\n";
            ofile << "endloop\n";
            ofile << "endfacet\n";
        }
    }
    ofile << "endsolid\n";
    ofile.close();

}

void TestBoxInterBasedSefIntersectionDetection(TriMesh& trimesh, const std::string& outfilename)
{
    auto start = std::chrono::system_clock::now();
    auto vecSelfInterPair = cglab::GetSelfIntersectionPairsBoxInterBased(trimesh);
    auto finish = std::chrono::system_clock::now();

    std::cout << "Box inter based get self inter pairs cost: " << std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count() << "ms" << std::endl;

    WriteOutInterPairs(outfilename, vecSelfInterPair, trimesh);
}

void TestAABspTreeBasedSelfIntersectionDetection(TriMesh& trimesh, const std::string& outfilename)
{
    auto start = std::chrono::system_clock::now();
    auto vecSelfInterPair = cglab::GetSelfIntersectionPairsAABspTreeBased(trimesh);
    auto finish = std::chrono::system_clock::now();

    std::cout << "AABspTree based get self inter pairs cost: " << std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count() << "ms" << std::endl;

    WriteOutInterPairs(outfilename, vecSelfInterPair, trimesh);
}

int main()
{
    // load trimesh and calculate distance
    TriMesh trimesh;
    //OpenMesh::IO::read_mesh(trimesh, std::string(MESHES_PATH) + "rabit.off");
    OpenMesh::IO::read_mesh(trimesh, "E:/STL/liver/RXF_20170608_133757-vein.obj");
    
    trimesh.request_face_normals();
    trimesh.request_vertex_normals();
    trimesh.update_normals();

    for (auto vh : trimesh.vertices())
    {
        auto normal = trimesh.normal(vh);
        auto& pt0 = trimesh.point(vh);
        pt0 += normal * 4;
    }
    OpenMesh::IO::write_mesh(trimesh, "D:/vein-offset.stl");

    TestBoxInterBasedSefIntersectionDetection(trimesh, "D:/vein-inter-pair-box-inter.stl");
    TestAABspTreeBasedSelfIntersectionDetection(trimesh, "D:/vein-inter-pair-bsp.stl");


    return 0;
}
