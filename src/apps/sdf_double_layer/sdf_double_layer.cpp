#include <iostream>
#include <fstream>

#include "openmesh_ext/openmesh_utils.h"
#include "config.h"

/*
* Wu, Y., Man, J., & Xie, Z. (2014). A double layer method for constructing signed distance fields from triangle meshes. Graphical Models, 76(4), 214–223. https://doi.org/10.1016/j.gmod.2014.04.011
*/
const unsigned int knmax = 16;

void SamplePoints(const TriMesh& trimesh, OpenMesh::FPropHandleT<double>& areaProperty, double epsilon, double pointOffsetSigma,
    std::vector<TriMesh::Point>& points, std::vector<TriMesh::Point>& pointsInner, std::vector<TriMesh::Point>& pointsExternal)
{
    int nface = trimesh.n_faces();

    for (int fidx = 0; fidx < nface; ++fidx)
    {
        unsigned int n = static_cast<unsigned int>(std::ceil(std::sqrt(trimesh.property(areaProperty, trimesh.face_handle(fidx)) / epsilon)));
        if (n > knmax)
        {
            n = knmax;
        }

        auto fviter = trimesh.cfv_iter(trimesh.face_handle(fidx));
        auto& ptA = trimesh.point(fviter++);
        auto& ptB = trimesh.point(fviter++);
        auto& ptC = trimesh.point(fviter++);

        TriMesh::Normal normaloffset = trimesh.normal(trimesh.face_handle(fidx)) * pointOffsetSigma;

        TriMesh::Point s1 = (ptB - ptA) / n;
        TriMesh::Point s2 = (ptC - ptB) / n;
        TriMesh::Point s3 = (ptC - ptA) / n;
        TriMesh::Point ps = 1.0 / 3.0 * (ptA + ptA + s1 + ptA + s3);
        points.push_back(ps);
        pointsInner.push_back(ps - normaloffset);
        pointsExternal.push_back(ps + normaloffset);
        unsigned int i = 2;
        while (i <= n)
        {
            ps = ps + s1;
            points.push_back(ps);
            pointsInner.push_back(ps - normaloffset);
            pointsExternal.push_back(ps + normaloffset);
            auto p = ps;
            unsigned int j = 2;
            while (j <= i)
            {
                p = p + s2;
                points.push_back(p);
                pointsInner.push_back(p - normaloffset);
                pointsExternal.push_back(p + normaloffset);
                j++;
            }
            i++;
        }
    }
}

int main()
{
    // load trimesh and calculate distance
    TriMesh trimesh;
    OpenMesh::IO::read_mesh(trimesh, std::string(MESHES_PATH) + "rabit.off");
    trimesh.request_face_normals();
    trimesh.update_normals();

    auto start = std::chrono::steady_clock::now();

    OpenMesh::FPropHandleT<double> faceAreaProperty;
    trimesh.add_property(faceAreaProperty);
    for (auto fh : trimesh.faces())
    {
        auto fviter = trimesh.fv_iter(fh);
        auto& pt0 = trimesh.point(fviter++);
        auto& pt1 = trimesh.point(fviter++);
        auto& pt2 = trimesh.point(fviter++);
        double area = 0.5 * std::abs(((pt0 - pt1).cross(pt0 - pt2)).norm());
        trimesh.property(faceAreaProperty, fh) = area;
    }

    TriMesh::Point low, high;
    GetTriMeshBoundingBox(trimesh, low, high);

    TriMesh::Point offset(0.1, 0.1, 0.1);
    low = low - offset;
    high = high + offset;

    int dim[] = { 512,512,512 };
    TriMesh::Point spacing = (high - low) / (TriMesh::Point(dim[0] - 1, dim[1] - 1, dim[2] - 1));

    // generate inner external points
    std::vector<TriMesh::Point> points;
    std::vector<TriMesh::Point> pointInner, pointExternal;
    points.reserve(trimesh.n_faces() * 2);
    pointInner.reserve(points.size());
    pointExternal.reserve(points.size());
    double areaEpsilon = std::max(std::max(spacing[0] * spacing[1], spacing[0] * spacing[2]), spacing[1]*spacing[2]);
    double pointOffsetSigma = std::max(std::max(spacing[0], spacing[1]), spacing[2]) * 2;
    SamplePoints(trimesh, faceAreaProperty, areaEpsilon, pointOffsetSigma, points, pointInner, pointExternal);
    auto end = std::chrono::steady_clock::now();

    std::cout << "sample points cost: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        << " ms" << std::endl;


    // compute internal external distance field
    int k = 4;
    int xySize = dim[0] * dim[1];
    int volumesize = xySize * dim[2];
    float* pInternalDS = new float[volumesize];
    memset(pInternalDS, 0x7f, volumesize * sizeof(float));
    //float* pExternalDS = new float[volumesize];

    float* xDisArray = new float[dim[0]];
    float* yDisArray = new float[dim[1]];
    float* zDisArray = new float[dim[2]];

    for (int i = 0; i < dim[0]; ++i)
    {
        xDisArray[i] = i * spacing[0] + low[0];
    }

    for (int i = 0; i < dim[1]; ++i)
    {
        yDisArray[i] = i * spacing[1] + low[1];
    }

    for (int i = 0; i < dim[2]; ++i)
    {
        zDisArray[i] = i * spacing[2] + low[2];
    }
    
    auto spacingInv = TriMesh::Point(1.0, 1.0, 1.) / spacing;
//#pragma omp parallel for num_threads(8)
    for (int i = 0; i < pointInner.size(); ++i)
    {
        auto offset = (pointInner[i] - low) * spacingInv;

        int x = (static_cast<int>(offset[0]));
        int y = (static_cast<int>(offset[1]));
        int z = (static_cast<int>(offset[2]));

        int xlow = std::max(x - k, 0);
        int ylow = std::max(y - k, 0);
        int zlow = std::max(z - k, 0);

        int xhigh = std::min(x + k, dim[0] - 1);
        int yhigh = std::min(y + k, dim[1] - 1);
        int zhigh = std::min(z + k, dim[2] - 1);

        for (int zidx = zlow; zidx <= zhigh; ++zidx)
        {
            float zDis = zDisArray[zidx] - pointInner[i][2];
            float zSquareDis = zDis * zDis;
            int zoffset = zidx * xySize;
            for (int yidx = ylow; yidx <= yhigh; ++yidx)
            {
                float yDis = yDisArray[yidx] - pointInner[i][1];
                float ySquareDis = yDis * yDis;
                int yoffset = yidx * dim[0];
                for (int xidx = xlow; xidx <= xhigh; ++xidx)
                {
                    float xDis = xDisArray[xidx] - pointInner[i][0];
                    float xSquareDis = xDis * xDis;
                    int idx = zoffset + yoffset + xidx;

                    pInternalDS[idx] = std::min(xSquareDis + ySquareDis + zSquareDis, pInternalDS[idx]);
                }
            }
        }
    }

    end = std::chrono::steady_clock::now();

    std::cout << "Generate volume cost: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        << " ms" << std::endl;

    std::ofstream ofile("D:/rabit.raw", std::ios::binary);
    ofile.write((char*)pInternalDS, volumesize * sizeof(float));
    ofile.close();

    delete[] xDisArray;
    delete[] yDisArray;
    delete[] zDisArray;
    delete[] pInternalDS;
    //delete[] pExternalDS;
    return 0;
}
