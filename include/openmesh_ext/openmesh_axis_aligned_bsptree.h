// ref: https://github.com/hahahuahai/axis-aligned-BSP/blob/master/bsp_tree.h

#ifndef OPENMESH_AXIS_ALIGNED_BSPTREE_H_
#define OPENMESH_AXIS_ALIGNED_BSPTREE_H_

#include <cstddef>
#include <vector>

#include "openmesh_typedef.h"
#include "openmesh_utils.h"

namespace cglab
{
namespace openmesh_ext
{

typedef int FaceIdx;

enum class Axis
{
    AXIS_X,
    AXIS_Y,
    AXIS_Z,
    AXIS_UNKNOW
};

struct AABspTreeNode
{
    std::vector<FaceIdx> faces;
    Axis axis = Axis::AXIS_UNKNOW;
    int depth = 0;
    float splitter = 0.0;
    TriMesh::Point minPoint;
    TriMesh::Point maxPoint;
    AABspTreeNode* leftChild = nullptr;
    AABspTreeNode* rightChild = nullptr;
};

class AABspTree
{
public:
    AABspTree() = default;
    ~AABspTree();

    void GenerateBspTree(TriMesh& triMesh, const TriMesh::Point& minPt, const TriMesh::Point& maxPt, int maxDepth, int nodeCapacity);
    const std::vector<AABspTreeNode *>& GetLeafCollection();
    void DeleteBspTree();

private:
    void DeleteTree_(AABspTreeNode* root);
    void SplitSpace_(AABspTreeNode* node, Axis axis, int depth);
    void CutFace_(FaceIdx face, Axis axis, float splitter, int& leftCount, int& rightCount, int& bothCount);

    std::vector<AABspTreeNode *> leafCollection_; // store all leaf nodes
    AABspTreeNode* root_ = nullptr;
    int maxDepth_ = 16;
    int nodeCapacity_ = 10;
    TriMesh triMesh_;
};

AABspTree::~AABspTree()
{
    if (root_ != nullptr)
    {
        DeleteBspTree();
    }
}

void AABspTree::GenerateBspTree(TriMesh& triMesh, const TriMesh::Point& minPt, const TriMesh::Point& maxPt, int maxDepth, int nodeCapacity)
{
    triMesh_ = triMesh;
    maxDepth_ = maxDepth;
    nodeCapacity_ = nodeCapacity;
    root_ = new AABspTreeNode();
    root_->depth = 1;
    root_->faces.reserve(triMesh.n_faces());
    for (auto f : triMesh.faces())
    {
        root_->faces.push_back(f.idx());
    }
    root_->leftChild = nullptr;
    root_->rightChild = nullptr;

    root_->minPoint = minPt;
    root_->maxPoint = maxPt;

    leafCollection_.reserve(triMesh.n_faces() / 6);
    SplitSpace_(root_, Axis::AXIS_UNKNOW, 1);
}

const std::vector<AABspTreeNode*>& AABspTree::GetLeafCollection()
{
    return leafCollection_;
}

void AABspTree::DeleteBspTree()
{
    DeleteTree_(root_);
    root_ = nullptr;
}

void AABspTree::DeleteTree_(AABspTreeNode* root)
{
    if (root != nullptr)
    {
        DeleteTree_(root->leftChild);
        DeleteTree_(root->rightChild);
        delete root;
        root = nullptr;
    }
}

void AABspTree::SplitSpace_(AABspTreeNode* node, Axis axis, int depth)
{
    if (node == nullptr)
    {
        return;
    }

    node->axis = axis;
    node->depth = depth;
    if (depth == maxDepth_)
    {
        leafCollection_.push_back(node);
        return;
    }

    if (node->faces.size() < nodeCapacity_)
    {
        leafCollection_.push_back(node);
        return;
    }

    node->leftChild = new AABspTreeNode();
    node->leftChild->faces.reserve(node->faces.size() / 2);
    node->rightChild = new AABspTreeNode();
    node->rightChild->faces.reserve(node->faces.size() / 2);

    node->leftChild->maxPoint = node->maxPoint;
    node->leftChild->minPoint = node->minPoint;
    node->rightChild->maxPoint = node->maxPoint;
    node->rightChild->minPoint = node->minPoint;

    float xLen = node->maxPoint[0] - node->minPoint[0];
    float yLen = node->maxPoint[1] - node->minPoint[1];
    float zLen = node->maxPoint[2] - node->minPoint[2];
    //设置该节点的划分点
    Axis mAxis = Axis::AXIS_X;
    if (yLen > xLen && yLen > zLen)
        mAxis = Axis::AXIS_Y;
    if (zLen > xLen && zLen > yLen)
        mAxis = Axis::AXIS_Z;
    switch (mAxis)
    {
    case Axis::AXIS_X:
        node->splitter = (node->maxPoint[0] + node->minPoint[0]) / 2;
        //取中间值	
        //修改子节点的包围盒值	
        node->leftChild->maxPoint[0] = node->splitter;
        node->rightChild->minPoint[0] = node->splitter;
        break;
    case Axis::AXIS_Y:
        node->splitter = (node->maxPoint[1] + node->minPoint[1]) / 2;
        node->leftChild->maxPoint[1] = node->splitter;
        node->rightChild->minPoint[1] = node->splitter;
        break;
    case Axis::AXIS_Z:
        node->splitter = (node->maxPoint[2] + node->minPoint[2]) / 2;
        node->leftChild->maxPoint[2] = node->splitter;
        node->rightChild->minPoint[2] = node->splitter;
        break;
    }

    for (int i = 0; i < node->faces.size(); ++i)
    {
        int leftCount = 0;
        int rightCount = 0;
        int bothCount = 0;
        //对每个面做切分到子节点	
        CutFace_(node->faces[i], mAxis, node->splitter, leftCount, rightCount, bothCount);
        if (leftCount || bothCount)
            node->leftChild->faces.push_back(node->faces[i]);
        if (rightCount || bothCount)
            node->rightChild->faces.push_back(node->faces[i]);
    }

    //关键的一步，分割完了之后，把当前层节点的面片清空	
    std::vector<FaceIdx>().swap(node->faces);
    
    //递归
    SplitSpace_(node->leftChild, mAxis, depth + 1);	
    SplitSpace_(node->rightChild, mAxis, depth + 1);
}

void AABspTree::CutFace_(FaceIdx face, Axis axis, float splitter, int& leftCount, int& rightCount, int& bothCount)
{
    auto fh = triMesh_.face_handle(face);
    auto fviter = triMesh_.fv_iter(fh);
    auto& pt0 = triMesh_.point(fviter++);
    auto& pt1 = triMesh_.point(fviter++);
    auto& pt2 = triMesh_.point(fviter++);
    
    float p[3]{ 0,0,0 };//记录某个方向的点分量
    switch (axis)
    {
    case Axis::AXIS_X:
        p[0] = pt0[0];
        p[1] = pt1[0];
        p[2] = pt2[0];
        break;
    case Axis::AXIS_Y:
        p[0] = pt0[1];
        p[1] = pt1[1];
        p[2] = pt2[1];
        break;
    case Axis::AXIS_Z:
        p[0] = pt0[2];
        p[1] = pt1[2];
        p[2] = pt2[2];
        break;
    }
    //对分隔到两边的顶点计数	
    for (int i = 0; i < 3; ++i)
    {
        if (p[i] < splitter)
            leftCount++;
        else if (p[i] > splitter)
            rightCount++;
        else if (fabs(p[i] - splitter) < 1e-6)
            bothCount++;
    }
}

} // namespace openmesh_ext
} // namespace cglab

#endif