#include <limits>
#include <set>
#include <vector>

#include "openmesh_ext/openmesh_typedef.h"

#ifdef USE_BOOST_RANDOM
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#endif

#include "Mathematics/IntrSegment3Triangle3.h"
#include "Mathematics/IntrTriangle3Triangle3.h"

namespace cglab
{
    namespace detail
    {
        const float inf = -std::numeric_limits<float>::max();
        const float sup = std::numeric_limits<float>::max();

        struct Box
        {
            typedef std::ptrdiff_t BoxID;
            float lo[3] = { sup, sup, sup };
            float hi[3] = { inf, inf, inf };

            int f_idx;
            BoxID ID() const
            {
                return (BoxID)(this);
            }
        };

        void UpdateBox(Box& box, const TriMesh::Point& pt)
        {
            for (int i = 0; i < 3; ++i)
            {
                if (pt[i] < box.lo[i])
                    box.lo[i] = pt[i];
                if (pt[i] > box.hi[i])
                    box.hi[i] = pt[i];
            }
        }

        struct BoxIntersectionCallback
        {
            void operator()(const Box* a, const Box* b)
            {
                if (a->ID() != b->ID())
                {
                    std::pair<int, int> facepair;
                    if (a->f_idx < b->f_idx)
                    {
                        facepair.first = a->f_idx;
                        facepair.second = b->f_idx;
                    }
                    else
                    {
                        facepair.first = b->f_idx;
                        facepair.second = a->f_idx;
                    }

                    auto& sharedSet = mapFaceIgnored[facepair.first];
                    if (sharedSet.find(facepair.second) == sharedSet.end())
                    {

                        vecIntersectionFaces.insert(facepair);
                    }
                }
            }

            std::set<std::pair<int, int>> vecIntersectionFaces;
            std::map<int, std::set<int>> mapFaceIgnored;
        };

        struct BoxPredicate
        {
            static auto CompareFunc(int dim)
            {
                return [=](const Box* a, const Box* b) -> bool {
                    return IsLoLessLo(a, b, dim);
                };
            }

            static auto SpanningFunc(double lo, double hi, int dim)
            {
                // 如果[lo, hi)在box的范围内，那么返回true
                return [=](const Box* a) -> bool {
                    return a->lo[dim] < lo && a->hi[dim] > hi;
                };
            }

            static auto LoLessFunc(double value, int dim)
            {
                return [=](const Box* a) -> bool {
                    return a->lo[dim] < value;
                };
            }

            static auto HiGreaterFunc(double value, int dim)
            {
                return [=](const Box* a) -> bool {
                    return HiGreater(a->hi[dim], value);
                };
            }

            static bool HiGreater(double hi, double val)
            {
                return hi >= val;
            }

            static bool IsLoLessLo(const Box* a, const Box* b, int dim)
            {
                return (a->lo[dim] < b->lo[dim]) ||
                    (a->lo[dim] == b->lo[dim] && ID(a) < ID(b));
            }

            static bool IsLoLessHi(const Box* a, const Box* b, int dim)
            {
                return HiGreater(b->hi[dim], a->lo[dim]);
            }

            static bool IsIntersect(const Box* a, const Box* b, int dim)
            {
                return IsLoLessHi(a, b, dim) && IsLoLessHi(b, a, dim);
            }

            static bool ContainsLoPoint(const Box* a, const Box* b, int dim)
            {
                return IsLoLessLo(a, b, dim) && IsLoLessHi(b, a, dim);
            }

            static size_t ID(const Box* a)
            {
                return reinterpret_cast<size_t>(a);
            }
        };

#ifdef USE_BOOST_RANDOM
        namespace mcsf_iterative_radon
        {
            template< class RandomAccessIter>
            RandomAccessIter Median_Of_Three(RandomAccessIter a, RandomAccessIter b, RandomAccessIter c, int dim)
            {
                if (BoxPredicate::IsLoLessLo(*a, *b, dim))
                {
                    if (BoxPredicate::IsLoLessLo(*b, *c, dim))
                    {
                        return b;
                    }
                    else if (BoxPredicate::IsLoLessLo(*a, *c, dim))
                    {
                        return c;
                    }
                    else
                    {
                        return a;
                    }
                }
                else if (BoxPredicate::IsLoLessLo(*a, *c, dim))
                {
                    return a;
                }
                else if (BoxPredicate::IsLoLessLo(*b, *c, dim))
                {
                    return c;
                }
                else
                {
                    return b;
                }
            }

            template< class RandomAccessIter>
            class Iterative_Radon
            {
                RandomAccessIter begin;
                std::ptrdiff_t size;
                int dim;

                boost::rand48 rng;
                boost::uniform_int<std::ptrdiff_t> dist;
                boost::variate_generator<boost::rand48&, boost::uniform_int<std::ptrdiff_t>> generator;

            public:
                Iterative_Radon(RandomAccessIter begin, RandomAccessIter end, int dim)
                    : begin(begin), size(end - begin), dim(dim)
                    , rng(), dist(0, size - 1), generator(rng, dist) {}

                RandomAccessIter operator () (int num_levels)
                {
                    if (num_levels < 0)
                    {
                        const std::ptrdiff_t d = generator();
                        return begin + d;
                    }

                    return Median_Of_Three(
                        (*this)(num_levels - 1),
                        (*this)(num_levels - 1),
                        (*this)(num_levels - 1), dim);
                }
            };
        }
#endif

        typedef std::function<void(const Box&, const Box&)> InterpCallback;

        template< class RandomAccessIter1, class RandomAccessIter2,
            class Callback>
            void OneWayScan(
                RandomAccessIter1 pBegin, RandomAccessIter1 pEnd,
                RandomAccessIter2 iBegin, RandomAccessIter2 iEnd,
                Callback& callback, int dim, bool in_order = true)
        {
            /// 对P和对I（根据I的low endpoint）进行排序；
            std::sort(pBegin, pEnd, BoxPredicate::CompareFunc(0));
            std::sort(iBegin, iEnd, BoxPredicate::CompareFunc(0));

            /// 遍历interval
            for (auto i = iBegin; i != iEnd; ++i)
            {
                /// 对P进行遍历，直到找到p1不小于第一个intervals的low endpoints点
                for (; pBegin != pEnd && BoxPredicate::IsLoLessLo((*pBegin), (*i), 0); ++pBegin) {}

                /// 继续对点进行遍历，直到找到一个点p2不小于第一个intervals的high endpoint点；
                /// 所有找到的点，属于第一个intervals；
                for (auto p = pBegin; p != pEnd && BoxPredicate::IsLoLessHi(*p, *i, 0); ++p)
                {
                    /// 如果第d为大于0，那么还需要直接判断其他维度是否相交，
                    /// 如果有某个维度不相交，那么i和p1就不相交
                    if ((std::ptrdiff_t)(*p) == (std::ptrdiff_t)(*i))
                        continue;
                    bool bIntersect = true;
                    for (int d = 1; d <= dim; ++d)
                    {
                        if (!BoxPredicate::IsIntersect(*p, *i, d))
                        {
                            bIntersect = false;
                            break;
                        }
                    }
                    if (bIntersect)
                    {
                        if (in_order)
                        {
                            callback(*p, *i);
                        }
                        else
                        {
                            callback(*i, *p);
                        }
                    }
                }
            }
        }

        template<class RandomAccessIter1, class RandomAccessIter2,
            class Callback>
            void ModifiedTwoWayScan(
                RandomAccessIter1 pBegin, RandomAccessIter1 pEnd,
                RandomAccessIter2 iBegin, RandomAccessIter2 iEnd,
                Callback& callback, int dim, bool in_order = true)
        {
            std::sort(pBegin, pEnd, BoxPredicate::CompareFunc(0));
            std::sort(iBegin, iEnd, BoxPredicate::CompareFunc(0));

            while (iBegin != iEnd && pBegin != pEnd)
            {
                // 此时i作为interval，p为points
                if (BoxPredicate::IsLoLessLo(*iBegin, *pBegin, 0))
                {
                    for (auto p = pBegin; p != pEnd && BoxPredicate::IsLoLessHi(*p, *iBegin, 0); ++p)
                    {
                        if ((std::ptrdiff_t)(*p) == (std::ptrdiff_t)(*iBegin))
                            continue;

                        bool bIntersect = true;
                        for (int d = 1; d <= dim; ++d)
                        {
                            if (!BoxPredicate::IsIntersect(*p, *iBegin, d))
                            {
                                bIntersect = false;
                                break;
                            }
                        }

                        if (bIntersect && BoxPredicate::ContainsLoPoint(*iBegin, *p, dim))
                        {
                            if (in_order)
                            {
                                callback(*p, *iBegin);
                            }
                            else
                            {
                                callback(*iBegin, *p);
                            }
                        }
                    }
                    ++iBegin;
                }
                else
                {
                    // 将p作为interval，i作为point
                    for (auto i = iBegin; i != iEnd && BoxPredicate::IsLoLessHi(*i, *pBegin, 0); ++i)
                    {
                        if ((std::ptrdiff_t)(*pBegin) == (std::ptrdiff_t)(*i))
                            continue;

                        bool bIntersect = true;
                        for (int d = 1; d <= dim; ++d)
                        {
                            if (!BoxPredicate::IsIntersect(*pBegin, *i, d))
                            {
                                bIntersect = false;
                                break;
                            }
                        }

                        if (bIntersect && BoxPredicate::ContainsLoPoint(*i, *pBegin, dim))
                        {
                            if (in_order)
                            {
                                callback(*pBegin, *i);
                            }
                            else
                            {
                                callback(*i, *pBegin);
                            }
                        }
                    }
                    ++pBegin;
                }
            }
        }

        template< class RandomAccessIter>
        RandomAccessIter SplitPoints(RandomAccessIter begin, RandomAccessIter end, int dim, double& mi)
        {
#ifdef USE_BOOST_RANDOM
            int levels = (int)(.91 * std::log(((double)std::distance(begin, end)) / 137.0) + 1);
            levels = (levels <= 0) ? 1 : levels;
            mcsf_iterative_radon::Iterative_Radon<RandomAccessIter> IR(begin, end, dim);
            RandomAccessIter it = IR(levels);
            mi = (*it)->lo[dim];
            return std::partition(begin, end, BoxPredicate::LoLessFunc(mi, dim));

#else
            std::ptrdiff_t delta = std::distance(begin, end);
            std::ptrdiff_t halfDelta = delta / 2;
            RandomAccessIter mid = begin + halfDelta;
            mi = (*mid)->lo[dim];
            return std::partition(begin, end, BoxPredicate::LoLessFunc(mi, dim));
#endif
        }

        template<class RandomAccessIter1, class RandomAccessIter2, class Callback>
        void SegmentTree(RandomAccessIter1 pBegin, RandomAccessIter1 pEnd,
            RandomAccessIter2 iBegin, RandomAccessIter2 iEnd,
            double lo, double hi, Callback& callback, std::ptrdiff_t cutoff, int dim, bool in_order)
        {

            if (pBegin == pEnd || iBegin == iEnd || lo >= hi)
            {
                return;
            }

            if (dim == 0)
            {
                //std::cout << "dim = 0. scanning ..." << std::endl;
                OneWayScan(pBegin, pEnd, iBegin, iEnd, callback, dim, in_order);
                return;
            }

            if (std::distance(pBegin, pEnd) < cutoff ||
                std::distance(iBegin, iEnd) < cutoff)
            {
                ModifiedTwoWayScan(pBegin, pEnd, iBegin, iEnd, callback, dim, in_order);
                return;
            }

            RandomAccessIter2 iSpanEnd = (lo == inf || hi == sup) ? iBegin :
                std::partition(iBegin, iEnd, BoxPredicate::SpanningFunc(lo, hi, dim));

            if (iBegin != iSpanEnd)
            {
                SegmentTree(pBegin, pEnd, iBegin, iSpanEnd, inf, sup, callback, cutoff, dim - 1, in_order);
                SegmentTree(iBegin, iSpanEnd, pBegin, pEnd, inf, sup, callback, cutoff, dim - 1, !in_order);
            }

            double mi(0);
            RandomAccessIter1 pMid = SplitPoints(pBegin, pEnd, dim, mi);

            if (pMid == pBegin || pMid == pEnd)
            {
                //std::cout << "unable to split points! performing modified two_way_scan ... " << std::endl;
                ModifiedTwoWayScan(pBegin, pEnd, iSpanEnd, iEnd, callback, dim, in_order);
                return;
            }

            RandomAccessIter2 iMid = std::partition(iSpanEnd, iEnd, BoxPredicate::LoLessFunc(mi, dim));
            //std::cout << "Processing ->left" << std::endl;
            SegmentTree(pBegin, pMid, iSpanEnd, iMid, lo, mi, callback, cutoff, dim, in_order);
            iMid = std::partition(iSpanEnd, iEnd, BoxPredicate::HiGreaterFunc(mi, dim));
            //std::cout << "Processing ->right" << std::endl;
            SegmentTree(pMid, pEnd, iSpanEnd, iMid, mi, hi, callback, cutoff, dim, in_order);
        }

        template<class RandomAccessIter1, class RandomAccessIter2, class Callback>
        void BoxIntersection(
            RandomAccessIter1 begin1, RandomAccessIter1 end1,
            RandomAccessIter2 begin2, RandomAccessIter2 end2,
            Callback& callback,
            std::ptrdiff_t cutoff = 10)
        {
            const int dim = 2; // start from 0

            SegmentTree(begin1, end1, begin2, end2, inf, sup, callback, cutoff, dim, true);
            // self intersection test, don't need below logic
            // SegmentTree(begin2, end2, begin1, end1, inf, sup, callback, cutoff, dim, false);
        }

        void ConvertMeshToBox(TriMesh& mesh, std::vector<Box>& vecBox, std::vector<Box*>& vecBoxPtr)
        {
            vecBox.reserve(mesh.n_faces());
            vecBoxPtr.reserve(mesh.n_faces());
            for (auto f : mesh.faces())
            {
                Box box;
                for (auto it_fv = mesh.fv_begin(f); it_fv != mesh.fv_end(f); it_fv++)
                {
                    auto point = mesh.point(*it_fv);
                    UpdateBox(box, point);
                }

                box.f_idx = f.idx();
                vecBox.push_back(box);
            }

            for (int i = 0; i < vecBox.size(); ++i)
            {
                vecBoxPtr.push_back(&vecBox[i]);
            }
        }
    } // namespace detail

    std::vector<std::pair<int, int>> GetSelfIntersectionPairsBoxInterBased(TriMesh& mesh)
    {
        /// 1. Box intersection test
        std::map<int, std::set<int>> mapFaceIgnored;
        std::map<std::pair<int, int>, int> mapFaceSharedVertexNotSharedEdge;

        // 先不考虑共顶点的三角形对
        // 共边三角形对
        for (auto f : mesh.faces())
        {
            std::set<int> setSharedEdges;
            for (auto ff_it = mesh.ff_begin(f); ff_it != mesh.ff_end(f); ff_it++)
            {
                setSharedEdges.insert((*ff_it).idx());
            }

            setSharedEdges.insert(f.idx());

            for (auto fv_it = mesh.fv_begin(f); fv_it != mesh.fv_end(f); fv_it++)
            {
                for (auto vf_it = mesh.vf_begin(fv_it); vf_it != mesh.vf_end(fv_it); vf_it++)
                {
                    auto aroundIdx = (*vf_it).idx();
                    if (f.idx() < aroundIdx)
                    {
                        mapFaceIgnored[f.idx()].insert(aroundIdx);
                        if (setSharedEdges.find(aroundIdx) == setSharedEdges.end())
                        {
                            std::pair<int, int> facepair{f.idx(), aroundIdx};
                            mapFaceSharedVertexNotSharedEdge[facepair] = (*fv_it).idx();
                        }
                    }
                }
            }
        }

        std::vector<detail::Box> vecMeshBoxA;
        std::vector<detail::Box*> vecMeshBoxPtrA;

        detail::ConvertMeshToBox(mesh, vecMeshBoxA, vecMeshBoxPtrA);

        detail::BoxIntersectionCallback callback;
        callback.mapFaceIgnored = mapFaceIgnored;
        std::vector<detail::Box*> vecMeshSelfBox(vecMeshBoxPtrA.begin(), vecMeshBoxPtrA.end());
        BoxIntersection(vecMeshBoxPtrA.begin(), vecMeshBoxPtrA.end(), vecMeshSelfBox.begin(), vecMeshSelfBox.end(), callback);

        std::set<std::pair<int, int>>& facePair = callback.vecIntersectionFaces;

        /// 2. Filter box intersection pairs by triangle intersection test
        std::vector<std::pair<int, int>> vecFilteredInterPair;
        vecFilteredInterPair.reserve(facePair.size() / 4);
        gte::TIQuery<float, gte::Triangle3<float>, gte::Triangle3<float>> queryTri3Tri3;
        for (auto iter = facePair.begin(); iter != facePair.end(); iter++)
        {
            gte::Triangle3<float> tri1, tri2;
            auto face1 = mesh.face_handle(iter->first);
            auto face2 = mesh.face_handle(iter->second);
            TriMesh::Point triPt[3];
            auto fviter = mesh.fv_iter(face1);
            triPt[0] = mesh.point(fviter++);
            triPt[1] = mesh.point(fviter++);
            triPt[2] = mesh.point(fviter++);

            tri1.v[0] = { triPt[0][0], triPt[0][1], triPt[0][2] };
            tri1.v[1] = { triPt[1][0], triPt[1][1], triPt[1][2] };
            tri1.v[2] = { triPt[2][0], triPt[2][1], triPt[2][2] };

            fviter = mesh.fv_iter(face2);
            triPt[0] = mesh.point(fviter++);
            triPt[1] = mesh.point(fviter++);
            triPt[2] = mesh.point(fviter++);

            tri2.v[0] = { triPt[0][0], triPt[0][1], triPt[0][2] };
            tri2.v[1] = { triPt[1][0], triPt[1][1], triPt[1][2] };
            tri2.v[2] = { triPt[2][0], triPt[2][1], triPt[2][2] };


            auto result = queryTri3Tri3(tri1, tri2);
            if (result.intersect)
            {
                vecFilteredInterPair.push_back(*iter);
            }
        }


        /// 3. Test the face pairs which sharing vertex and not sharing edge

        // Filter from box inter pair

        // 从共顶点三角形对中筛选出相交三角形对
        gte::TIQuery<float, gte::Segment3<float>, gte::Triangle3<float>> query;
        for (auto iter = mapFaceSharedVertexNotSharedEdge.begin(); iter != mapFaceSharedVertexNotSharedEdge.end(); ++iter)
        {
            auto interPair = iter->first;
            auto sharedVertexIdx = iter->second;

            auto vh = mesh.vertex_handle(sharedVertexIdx);
            auto f1h = mesh.face_handle(interPair.first);
            auto f2h = mesh.face_handle(interPair.second);
            auto f1hh = mesh.halfedge_handle(f1h);
            auto f2hh = mesh.halfedge_handle(f2h);

            while (mesh.to_vertex_handle(mesh.next_halfedge_handle(f1hh)) != vh)
            {
                f1hh = mesh.next_halfedge_handle(f1hh);
            }

            while (mesh.to_vertex_handle(mesh.next_halfedge_handle(f2hh)) != vh)
            {
                f2hh = mesh.next_halfedge_handle(f2hh);
            }

            // 计算顶点到两条对边的距离
            auto f1_edge_pt0 = mesh.point(mesh.from_vertex_handle(f1hh));
            auto f1_edge_pt1 = mesh.point(mesh.to_vertex_handle(f1hh));
            auto f2_edge_pt0 = mesh.point(mesh.from_vertex_handle(f2hh));
            auto f2_edge_pt1 = mesh.point(mesh.to_vertex_handle(f2hh));
            auto pt2 = mesh.point(vh);

            bool bInter = false;
            {
                // 计算f1的边f1 pt0 pt1是否和三角形f2相交
                gte::Segment3<float> line;
                line.p[0] = { f1_edge_pt0[0], f1_edge_pt0[1], f1_edge_pt0[2] };
                line.p[1] = { f1_edge_pt1[0], f1_edge_pt1[1], f1_edge_pt1[2] };

                gte::Triangle3<float> tri;
                tri.v[0] = { f2_edge_pt0[0], f2_edge_pt0[1], f2_edge_pt0[2] };
                tri.v[1] = { f2_edge_pt1[0], f2_edge_pt1[1], f2_edge_pt1[2] };
                tri.v[2] = { pt2[0], pt2[1], pt2[2] };
                
                auto result = query(line, tri);
                bInter = bInter || result.intersect;
            }
            {
                // 计算f2的边f2 pt0 pt1是否和三角形f1相交
                gte::Segment3<float> line;
                line.p[0] = { f2_edge_pt0[0], f2_edge_pt0[1], f2_edge_pt0[2] };
                line.p[1] = { f2_edge_pt1[0], f2_edge_pt1[1], f2_edge_pt1[2] };

                gte::Triangle3<float> tri;
                tri.v[0] = { f1_edge_pt0[0], f1_edge_pt0[1], f1_edge_pt0[2] };
                tri.v[1] = { f1_edge_pt1[0], f1_edge_pt1[1], f1_edge_pt1[2] };
                tri.v[2] = { pt2[0], pt2[1], pt2[2] };

                auto result = query(line, tri);
                bInter = bInter || result.intersect;
            }

            if (bInter)
            {
                vecFilteredInterPair.push_back(interPair);
            }
        }

        return vecFilteredInterPair;
    }

} // namespace cglab