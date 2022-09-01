// ref: https://github.com/brandonpelfrey/SimpleOctree
// https://www.cnblogs.com/Glucklichste/p/11505743.html

#ifndef OPENMESH_Octree_H
#define OPENMESH_Octree_H

#include <cstddef>
#include <vector>

#include "openmesh_typedef.h"
#include "openmesh_utils.h"

/**!
 * Octree data struct used for openmesh datastruct, the vertices of trimesh is stored in every oct-grid
 */
class OpenMeshOctree {
    /*
        Children follow a predictable pattern to make accesses simple.
        Here, - means less than 'origin' in that dimension, + means greater than.
        child:	0 1 2 3 4 5 6 7
        x:      - - - - + + + +
        y:      - - + + - - + +
        z:      - + - + - + - +
     */
public:
    OpenMeshOctree(const TriMesh& triMesh, const TriMesh::Point& origin, const TriMesh::Point& halfDimension, unsigned int currentDepth)
    : triMesh_(triMesh), origin_(origin), halfDimension_(halfDimension), currentDepth_(currentDepth)
    {
        // Initially, there are no children
        for(int i=0; i<8; ++i) 
            children_[i] = nullptr;
    }

    OpenMeshOctree(const OpenMeshOctree& copy) = delete;

    ~OpenMeshOctree() 
    {
        // Recursively destroy octants
        for (int i = 0; i < 8; ++i)
        {
            if (children_[i])
                delete children_[i];
        }
    }

    void Insert(const TriMesh::VertexHandle& vh) {
        if (currentDepth_ == 0)  // Reach the max depth
        {
            data_.insert(vh);
            return;
        }
        else if (!hasChild_ && data_.empty()) // the current node has no data and has no child
        {
            data_.insert(vh);
            return;
        }
        else if (!hasChild_ && !data_.empty()) // the current node has no child and has data
        {
            // We're at a leaf, but there's already something here
            // We will split this node so that it has 8 child octants
            // and then insert the old data that was here, along with 
            // this new data point
            std::set<TriMesh::VertexHandle> oldData;
            oldData.swap(data_);

            assert(oldData.size() == 1);

            // only create used child
            int idxChild0 = GetOctantContainingPoint(*oldData.begin());
            int idxChild1 = GetOctantContainingPoint(vh);
            int idxChild[] = { idxChild0, idxChild1 };
            int idxChildCount = idxChild0 == idxChild1 ? 1 : 2;
            for (int i = 0; i < idxChildCount; ++i)
            {
                if (children_[idxChild[i]] == nullptr)
                {
                    TriMesh::Point newOrigin = origin_;
                    newOrigin[0] += halfDimension_[0] * (idxChild[i] & 4 ? .5f : -.5f);
                    newOrigin[1] += halfDimension_[1] * (idxChild[i] & 2 ? .5f : -.5f);
                    newOrigin[2] += halfDimension_[2] * (idxChild[i] & 1 ? .5f : -.5f);
                    children_[idxChild[i]] = new OpenMeshOctree(triMesh_, newOrigin, halfDimension_ * .5f, currentDepth_ - 1);
                }
            }

            // Re-insert the old point, and insert this new point
            // (We wouldn't need to insert from the root, because we already
            // know it's guaranteed to be in this section of the tree)
            children_[idxChild0]->Insert(*oldData.begin());
            children_[idxChild1]->Insert(vh);
            hasChild_ = true;
        }
        else
        {
            // We are at an interior node. Insert recursively into the 
            // appropriate child octant
            int octant = GetOctantContainingPoint(vh);
            if (children_[octant] == nullptr)
            {
                TriMesh::Point newOrigin = origin_;
                newOrigin[0] += halfDimension_[0] * (octant & 4 ? .5f : -.5f);
                newOrigin[1] += halfDimension_[1] * (octant & 2 ? .5f : -.5f);
                newOrigin[2] += halfDimension_[2] * (octant & 1 ? .5f : -.5f);
                children_[octant] = new OpenMeshOctree(triMesh_, newOrigin, halfDimension_ * .5f, currentDepth_ - 1);
            }
            children_[octant]->Insert(vh);
        }
    }

    // The PointType must support
    // pt[0], pt[1], pt[2] to access the coordinate value
    template <typename PointType>
    void GetNodeContainPoint(const PointType& pt, OpenMeshOctree*& result)
    {
        // Compute the min/max corners of this child octant
        auto cmax = origin_ + halfDimension_;
        auto cmin = origin_ - halfDimension_;

        // If the query rectangle is outside the child's bounding box, 
        // then continue
        if (cmax[0] < pt[0] || cmax[1] < pt[1] || cmax[2] < pt[2]) return;
        if (cmin[0] > pt[0] || cmin[1] > pt[1] || cmin[2] > pt[2]) return;

        result = this;

        bool isLeafNode = !hasChild_ || currentDepth_ == 0;
        if (!isLeafNode)
        {
            for (int i = 0; i < 8; ++i) {
                if (children_[i] == nullptr) continue;

                // At this point, we've determined that this child is containing the pt 
                children_[i]->GetNodeContainPoint(pt, result);
            }
        }
    }

    void GetDataInCurrentNode(std::set<TriMesh::VertexHandle>& data)
    {
        bool isLeafNode = !hasChild_ || currentDepth_ == 0;

        if (isLeafNode)
        {
            data.insert(data_.begin(), data_.end());
        }
        else
        {
            for (int i = 0; i < 8; ++i)
            {
                if (children_[i] == nullptr) continue;
                children_[i]->GetDataInCurrentNode(data);
            }
        }
    }

private:
    // Determine which octant of the tree would contain 'point'
    int GetOctantContainingPoint(const TriMesh::VertexHandle& vh) const
    {
        const TriMesh::Point& point = triMesh_.point(vh);
        int oct = 0;
        if (point[0] >= origin_[0]) oct |= 4;
        if (point[1] >= origin_[1]) oct |= 2;
        if (point[2] >= origin_[2]) oct |= 1;
        return oct;
    }

    const TriMesh& triMesh_;

    // Physical position/size. This implicitly defines the bounding 
    // box of this node
    TriMesh::Point origin_;         //! The physical center of this node
    TriMesh::Point halfDimension_;  //! Half the width/height/depth of this node

    // The tree has up to eight children and can additionally store
    // a point, though in many applications only, the leaves will store data.
    OpenMeshOctree* children_[8]; //! Pointers to child octants
    std::set<TriMesh::VertexHandle> data_;   //! Data point to be stored at a node

    unsigned int currentDepth_ = 0; // The depth of root is max, the depth of leaf node is 0 or all child is null
    bool hasChild_ = false;
};

/// \brief Create octree data struct from openmesh trimesh 
/// \param[in] trimesh input triangle mesh
/// \param[in] maxDepth max depth of octree, the depth is from 16(root) to 0
OpenMeshOctree* CreateOpenMeshOctree(const TriMesh& trimesh, unsigned int maxDepth = 16)
{
    TriMesh::Point low, high;
    GetTriMeshBoundingBox(trimesh, low, high);

    TriMesh::Point center = (low + high) * 0.5;
    TriMesh::Point halfDimension = (high - low) * 0.5;

    OpenMeshOctree* root = new OpenMeshOctree(trimesh, center, halfDimension, maxDepth);
    for (auto vh : trimesh.vertices())
    {
        root->Insert(vh);
    }

    return root;
}

OpenMeshOctree* CreateOpenMeshOctree(const TriMesh& trimesh, const TriMesh::Point& low, const TriMesh::Point& high, unsigned int maxDepth = 16)
{
    TriMesh::Point center = (low + high) * 0.5;
    TriMesh::Point halfDimension = (high - low) * 0.5;

    OpenMeshOctree* root = new OpenMeshOctree(trimesh, center, halfDimension, maxDepth);
    for (auto vh : trimesh.vertices())
    {
        root->Insert(vh);
    }

    return root;
}

#endif