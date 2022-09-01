#include <osg/Texture1D>
#include <osg/Geode>
#include <osgViewer/Viewer>

#include "openmesh_ext/openmesh_octree.h"
#include "openmesh_ext/openmesh_to_osg_mesh.h"
#include "Mathematics/DistPointTriangle.h"

#include "common/colormap_define.h"
#include "common/osg_common.h"

double GetPointToFaceDistance(TriMesh& trimesh1, TriMesh::VertexHandle vh, TriMesh& trimesh2, TriMesh::FaceHandle fh)
{
    auto pt = trimesh1.point(vh);
    gte::DCPPoint3Triangle3<double> distanceQuery;
    gte::Vector<3, double> point = { pt[0], pt[1], pt[2] };

    TriMesh::Point triPt[3];
    auto fviter = trimesh2.fv_iter(fh);
    triPt[0] = trimesh2.point(fviter++);
    triPt[1] = trimesh2.point(fviter++);
    triPt[2] = trimesh2.point(fviter++);

    gte::Triangle<3, double> triangle;
    triangle.v[0] = {triPt[0][0], triPt[0][1], triPt[0][2]};
    triangle.v[1] = {triPt[1][0], triPt[1][1], triPt[1][2]};
    triangle.v[2] = {triPt[2][0], triPt[2][1], triPt[2][2]};

    return distanceQuery(point, triangle).distance;
}

osg::Camera* CreateHUDContent(double left, double right, double bottom, double top, double minV, double maxV)
{
    osg::ref_ptr<osg::Camera> hudcamera = CreateHUDCamera(left, right, bottom, top);

    std::stringstream ss1;
    ss1.precision(2);
    ss1.setf(std::ios::fixed);
    ss1 << "minValue " << minV;

    std::string sMin = ss1.str();
    ss1.clear();
    ss1.str("");

    ss1 << "maxValue " << maxV;
    std::string sMax = ss1.str();

    osg::ref_ptr<osgText::Text> text1 = CreateText(osg::Vec3(600, 50, 0), sMin, 20);
    osg::ref_ptr<osgText::Text> text2 = CreateText(osg::Vec3(600, 550, 0), sMax, 20);

    osg::ref_ptr<osg::Geode> geode = new osg::Geode;
    geode->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
    geode->addDrawable(text1);
    geode->addDrawable(text2);
    geode->addDrawable(CreateColorMap(right - 50, right, bottom+5, top-5));
    hudcamera->addChild(geode);
    return hudcamera.release();
}

int main()
{
    // load trimesh and calculate distance
    TriMesh trimeshbase, trimesh;
    OpenMesh::IO::read_mesh(trimeshbase, "E:/STL/test/liver.obj");
    OpenMesh::IO::read_mesh(trimesh, "E:/STL/test/liver-out-0.5.stl");

    TriMesh::Point low, high;
    GetTriMeshBoundingBox(trimesh, low, high);

    TriMesh::Point lowbase, highbase;
    GetTriMeshBoundingBox(trimeshbase, lowbase, highbase);

    low[0] = low[0] < lowbase[0] ? low[0] : lowbase[0];
    low[1] = low[1] < lowbase[1] ? low[1] : lowbase[1];
    low[2] = low[2] < lowbase[2] ? low[2] : lowbase[2];

    high[0] = high[0] > highbase[0] ? high[0] : highbase[0];
    high[1] = high[1] > highbase[1] ? high[1] : highbase[1];
    high[2] = high[2] > highbase[2] ? high[2] : highbase[2];

    TriMesh::Point offset(0.1, 0.1, 0.1);

    low = low - offset;
    high = high + offset;

    OpenMeshOctree* octree = CreateOpenMeshOctree(trimeshbase, low, high);

    OpenMesh::VPropHandleT<double>  vertexDistanceProperty;

    trimesh.add_property(vertexDistanceProperty);
    
    double totalminDis = 10e30;
    double totalmaxDis = -100;

    // find the nearst point
    OpenMeshOctree* node = nullptr;
    for (auto vh : trimesh.vertices())
    {
        octree->GetNodeContainPoint(trimesh.point(vh), node);
        if (node != nullptr)
        {
            std::set<TriMesh::VertexHandle> vertices;
            node->GetDataInCurrentNode(vertices);
            std::set<TriMesh::FaceHandle> faces;
            for (auto iterV = vertices.begin(); iterV != vertices.end(); ++iterV)
            {
                for (auto iterF = trimeshbase.vf_begin(*iterV); iterF != trimeshbase.vf_end(*iterV); ++iterF)
                {
                    faces.insert(*iterF);
                }
            }

            assert(!vertices.empty());
            
            double minDistance = 10e30;
            for (auto fh : faces)
            {
                double ptToFaceDis = GetPointToFaceDistance(trimesh, vh, trimeshbase, fh);
                if (ptToFaceDis < minDistance)
                {
                    minDistance = ptToFaceDis;
                }
            }

            trimesh.property(vertexDistanceProperty, vh) = minDistance;
            if (totalmaxDis < minDistance) totalmaxDis = minDistance;
            if (totalminDis > minDistance) totalminDis = minDistance;
        }
    }

    // construct osg geometry
    osg::ref_ptr<osg::Geometry> geom = new osg::Geometry;
    TriMeshToGeometry(trimesh, geom);

    // Calculate color on cpu
    double red, green, blue;
    osg::ref_ptr<osg::Vec3Array> colors = new osg::Vec3Array(trimesh.n_vertices());
    for (int i = 0; i < trimesh.n_vertices(); ++i)
    {
        double& dis = trimesh.property(vertexDistanceProperty, trimesh.vertex_handle(i));
        InterplateColor((dis - totalminDis) / (totalmaxDis - totalminDis), red, green, blue);

        (*colors)[i].set(red, green, blue);
    }
    geom->setColorArray(colors.get());
    geom->setColorBinding(osg::Geometry::BIND_PER_VERTEX);

    osg::ref_ptr<osg::Geode> geode = new osg::Geode;
    geode->addDrawable(geom);

    osg::ref_ptr<osg::Group> root = new osg::Group;
    root->addChild(geode);

    trimesh.remove_property(vertexDistanceProperty);
    delete octree;

    root->addChild(CreateHUDContent(0, 800, 0, 600, totalminDis, 0.5));

    // transform by camera
    osgViewer::Viewer viewer;
    viewer.setSceneData(root);
    viewer.setUpViewInWindow(100, 100, 800, 600);
    return viewer.run();
}
