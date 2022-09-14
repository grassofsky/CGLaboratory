#include <iostream>
#include <fstream>

#include "openmesh_ext/openmesh_utils.h"
#include "config.h"

const unsigned int knmax = 16;

int main()
{
    // load trimesh and calculate distance
    TriMesh trimesh;
    OpenMesh::IO::read_mesh(trimesh, std::string(MESHES_PATH) + "rabit.off");
    trimesh.request_face_normals();
    trimesh.update_normals();

    return 0;
}
