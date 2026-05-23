#ifndef ASN_MESH_TET_MESH_HPP
#define ASN_MESH_TET_MESH_HPP

#include <iomanip>
#include <iostream>
#include <fstream>

#include "TetMesh.h"

#include "Utils.hpp"
#include "Tree.hpp"

namespace Asn::Mesh
{
    inline const List<Pair<Array<Real, 4>, Array<Int, 2>>> TET_TWICE_REFINE_TABLE = {
                {{0.75, 0.25, 0.00, 0.00}, {0, 4}},
                {{0.75, 0.00, 0.25, 0.00}, {0, 5}},
                {{0.75, 0.00, 0.00, 0.25}, {0, 6}},

                {{0.25, 0.75, 0.00, 0.00}, {1, 4}},
                {{0.00, 0.75, 0.25, 0.00}, {1, 7}},
                {{0.00, 0.75, 0.00, 0.25}, {1, 8}},

                {{0.25, 0.00, 0.75, 0.00}, {2, 5}},
                {{0.00, 0.25, 0.75, 0.00}, {2, 7}},
                {{0.00, 0.00, 0.75, 0.25}, {2, 9}},

                {{0.25, 0.00, 0.00, 0.75}, {3, 6}},
                {{0.00, 0.25, 0.00, 0.75}, {3, 8}},
                {{0.00, 0.00, 0.25, 0.75}, {3, 9}},

                {{0.50, 0.25, 0.25, 0.00}, {4, 5}},
                {{0.25, 0.50, 0.25, 0.00}, {4, 7}},
                {{0.25, 0.25, 0.50, 0.00}, {5, 7}},

                {{0.50, 0.25, 0.00, 0.25}, {4, 6}},
                {{0.25, 0.50, 0.00, 0.25}, {4, 8}},
                {{0.25, 0.25, 0.00, 0.50}, {6, 8}},

                {{0.50, 0.00, 0.25, 0.25}, {5, 6}},
                {{0.25, 0.00, 0.50, 0.25}, {5, 9}},
                {{0.25, 0.00, 0.25, 0.50}, {6, 9}},

                {{0.00, 0.50, 0.25, 0.25}, {7, 8}},
                {{0.00, 0.25, 0.50, 0.25}, {7, 9}},
                {{0.00, 0.25, 0.25, 0.50}, {8, 9}}
            };

    template <Real TOL>
    inline Order TetMesh<TOL>::cmp_vert_idx(const Int &i, const Int &j) const
    {
        if (CMP_POINT(vert[i].coord, vert[j].coord))
        {
            return Order::less;
        }
        if (CMP_POINT(vert[j].coord, vert[i].coord))
        {
            return Order::greater;
        }
        return Order::equal;
    }

    template <Real TOL>
    template <Int N>
    inline void TetMesh<TOL>::sort_by_pos(Array<Int, N> &list) const
    {
        for (Int i = 0; i < N; i++)
        {
            Int idx = i;
            for (Int j = i + 1; j < N; j++)
            {
                const Order order = cmp_vert_idx(list[idx], list[j]);
                assert(order != Order::equal);
                if (order == Order::greater)
                {
                    idx = j;
                }
            }
            std::swap(list[i], list[idx]);
        }
        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::load_mesh(const String &filename)
    {
        clear();

        std::ifstream fp(filename, std::ios::in);
        if (!fp.is_open())
        {
            ASN_ERROR("Can't open input file.");
        }

        Int num;

        if (!seek_keyword(fp, "Vertices"))
        {
            ASN_ERROR("Can't find keyword \"Vertices\".");
        }
        fp >> num;
        if (fp.fail() || num <= 0)
        {
            ASN_ERROR("Fail reading file.");
        }
        for (Int i = 0; i < num; i++)
        {
            Real x, y, z;
            Int l;
            fp >> x >> y >> z >> l;
            if (fp.fail())
            {
                ASN_ERROR("Fail reading file.");
            }
            Int vid = add_vert(Point<SPACE_DIM>({x, y, z}));
            vert[vid].label = Unknown;
        }

        if (!seek_keyword(fp, "Tetrahedra"))
        {
            ASN_ERROR("Can't find keyword \"Tetrahedra\".");
        }
        fp >> num;
        if (fp.fail() || num < 0)
        {
            ASN_ERROR("Fail reading file.");
        }
        for (Int i = 0; i < num; i++)
        {
            Int p0, p1, p2, p3, l;
            fp >> p0 >> p1 >> p2 >> p3 >> l;
            if (fp.fail())
            {
                ASN_ERROR("Fail reading file.");
            }
            if (p0 < 1 || p1 < 1 || p2 < 1 || p3 < 1)
            {
                ASN_ERROR("Fail reading file.");
            }
            Int pid = add_poly({p0 - 1, p1 - 1, p2 - 1, p3 - 1}, nullptr);
            poly[pid].label = Unknown;
        }
        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::save_mesh(const String &filename) const
    {
        std::ofstream fp(filename, std::ios::out);
        if (!fp.is_open())
        {
            ASN_ERROR("Can't open output file.");
        }
        fp << "MeshVersionFormatted 1" << std::endl;
        fp << "Dimension 3" << std::endl;
        fp << std::endl;
        Int num = vert.size();
        if (num > 0)
        {
            fp << "Vertices" << std::endl;
            fp << num << std::endl;
            for (Int vid = 0; vid < num; vid++)
            {
                const Vert &v = vert[vid];
                fp << std::setprecision(15) << v.coord(0) << " ";
                fp << std::setprecision(15) << v.coord(1) << " ";
                fp << std::setprecision(15) << v.coord(2) << " ";
                fp << v.label << std::endl;
            }
        }
        fp << std::endl;
        num = poly.size();
        if (num > 0)
        {
            fp << "Tetrahedra" << std::endl;
            fp << num << std::endl;
            for (Int pid = 0; pid < num; pid++)
            {
                const Poly &p = poly[pid];
                for (const Int id : p.vert)
                {
                    fp << (id + 1) << " ";
                }
                fp << p.label << std::endl;
            }
        }
        fp << std::endl;
        fp << "End" << std::endl;
        fp.close();
        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::save_obj(const String &filename) const
    {
        std::ofstream fp(filename, std::ios::out);
        if (!fp.is_open())
        {
            ASN_ERROR("Can't open output file.");
        }
        Int num = vert.size();
        for (Int vid = 0; vid < num; vid++)
        {
            const Vert &v = vert[vid];
            fp << "v ";
            fp << std::setprecision(15) << v.coord(0) << " ";
            fp << std::setprecision(15) << v.coord(1) << " ";
            fp << std::setprecision(15) << v.coord(2) << " ";
            fp << std::endl;
        }
        fp << std::endl;
        num = face.size();
        for (Int fid = 0; fid < num; fid++)
        {
            const Face &f = face[fid];
            fp << "f";
            for (const Int vid : f.vert)
            {
                fp << " " << (vid + 1);
            }
            fp << std::endl;
        }
        fp.close();
        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::clear()
    {
        dangling_face.clear();
        dangling_edge.clear();
        dangling_vert.clear();

        vert_index.clear();

        vert.clear();
        edge.clear();
        face.clear();
        poly.clear();

        tree.clear();

        red_refine_poly.clear();
        green_refine_poly.clear();
        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::switch_vert(const Int &vid1, const Int &vid2)
    {
        if (vid1 == vid2)
        {
            return;
        }

        assert(0 <= vid1 && vid1 < static_cast<Int>(vert.size()));
        assert(0 <= vid2 && vid2 < static_cast<Int>(vert.size()));

        std::swap(vert[vid1], vert[vid2]);
        std::swap(vert_index[vert[vid1].coord], vert_index[vert[vid2].coord]);

        HashSet<Int> verts_to_update;
        verts_to_update.insert(vert[vid1].vert.begin(), vert[vid1].vert.end());
        verts_to_update.insert(vert[vid2].vert.begin(), vert[vid2].vert.end());

        HashSet<Int> edges_to_update;
        edges_to_update.insert(vert[vid1].edge.begin(), vert[vid1].edge.end());
        edges_to_update.insert(vert[vid2].edge.begin(), vert[vid2].edge.end());

        HashSet<Int> faces_to_update;
        faces_to_update.insert(vert[vid1].face.begin(), vert[vid1].face.end());
        faces_to_update.insert(vert[vid2].face.begin(), vert[vid2].face.end());

        HashSet<Int> polys_to_update;
        polys_to_update.insert(vert[vid1].poly.begin(), vert[vid1].poly.end());
        polys_to_update.insert(vert[vid2].poly.begin(), vert[vid2].poly.end());

        for (Int i : verts_to_update)
        {
            assert(0 <= i && i < static_cast<Int>(vert.size()));
            switch_data(vert[i].vert, vid1, vid2);
        }

        for (Int i : edges_to_update)
        {
            assert(0 <= i && i < static_cast<Int>(edge.size()));
            switch_data_array<VERTS_PER_EDGE>(edge[i].vert, vid1, vid2);
        }

        for (Int i : faces_to_update)
        {
            assert(0 <= i && i < static_cast<Int>(face.size()));
            switch_data_array<VERTS_PER_FACE>(face[i].vert, vid1, vid2);
        }

        for (Int i : polys_to_update)
        {
            assert(0 <= i && i < static_cast<Int>(poly.size()));
            switch_data_array<VERTS_PER_POLY>(poly[i].vert, vid1, vid2);
        }
        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::switch_edge(const Int &eid1, const Int &eid2)
    {
        if (eid1 == eid2)
        {
            return;
        }

        assert(0 <= eid1 && eid1 < static_cast<Int>(edge.size()));
        assert(0 <= eid2 && eid2 < static_cast<Int>(edge.size()));

        std::swap(edge[eid1], edge[eid2]);

        HashSet<Int> verts_to_update;
        verts_to_update.insert(edge[eid1].vert.begin(), edge[eid1].vert.end());
        verts_to_update.insert(edge[eid2].vert.begin(), edge[eid2].vert.end());

        HashSet<Int> edges_to_update;
        edges_to_update.insert(edge[eid1].edge.begin(), edge[eid1].edge.end());
        edges_to_update.insert(edge[eid2].edge.begin(), edge[eid2].edge.end());

        HashSet<Int> faces_to_update;
        faces_to_update.insert(edge[eid1].face.begin(), edge[eid1].face.end());
        faces_to_update.insert(edge[eid2].face.begin(), edge[eid2].face.end());

        HashSet<Int> polys_to_update;
        polys_to_update.insert(edge[eid1].poly.begin(), edge[eid1].poly.end());
        polys_to_update.insert(edge[eid2].poly.begin(), edge[eid2].poly.end());

        for (Int i : verts_to_update)
        {
            assert(0 <= i && i < static_cast<Int>(vert.size()));
            switch_data(vert[i].edge, eid1, eid2);
        }

        for (Int i : edges_to_update)
        {
            assert(0 <= i && i < static_cast<Int>(edge.size()));
            switch_data(edge[i].edge, eid1, eid2);
        }

        for (Int i : faces_to_update)
        {
            assert(0 <= i && i < static_cast<Int>(face.size()));
            switch_data(face[i].edge, eid1, eid2);
        }

        for (Int i : polys_to_update)
        {
            assert(0 <= i && i < static_cast<Int>(poly.size()));
            switch_data(poly[i].edge, eid1, eid2);
        }

        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::switch_face(const Int &fid1, const Int &fid2)
    {
        if (fid1 == fid2)
        {
            return;
        }

        assert(0 <= fid1 && fid1 < static_cast<Int>(face.size()));
        assert(0 <= fid2 && fid2 < static_cast<Int>(face.size()));

        std::swap(face[fid1], face[fid2]);

        HashSet<Int> verts_to_update;
        verts_to_update.insert(face[fid1].vert.begin(), face[fid1].vert.end());
        verts_to_update.insert(face[fid2].vert.begin(), face[fid2].vert.end());

        HashSet<Int> edges_to_update;
        edges_to_update.insert(face[fid1].edge.begin(), face[fid1].edge.end());
        edges_to_update.insert(face[fid2].edge.begin(), face[fid2].edge.end());

        HashSet<Int> faces_to_update;
        faces_to_update.insert(face[fid1].face.begin(), face[fid1].face.end());
        faces_to_update.insert(face[fid2].face.begin(), face[fid2].face.end());

        HashSet<Int> polys_to_update;
        polys_to_update.insert(face[fid1].poly.begin(), face[fid1].poly.end());
        polys_to_update.insert(face[fid2].poly.begin(), face[fid2].poly.end());

        for (Int i : verts_to_update)
        {
            assert(0 <= i && i < static_cast<Int>(vert.size()));
            switch_data(vert[i].face, fid1, fid2);
        }

        for (Int i : edges_to_update)
        {
            assert(0 <= i && i < static_cast<Int>(edge.size()));
            switch_data(edge[i].face, fid1, fid2);
        }

        for (Int i : faces_to_update)
        {
            assert(0 <= i && i < static_cast<Int>(face.size()));
            switch_data(face[i].face, fid1, fid2);
        }

        for (Int i : polys_to_update)
        {
            assert(0 <= i && i < static_cast<Int>(poly.size()));
            switch_data(poly[i].face, fid1, fid2);
        }

        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::switch_poly(const Int &pid1, const Int &pid2)
    {
        if (pid1 == pid2)
        {
            return;
        }

        assert(0 <= pid1 && pid1 < static_cast<Int>(poly.size()));
        assert(0 <= pid2 && pid2 < static_cast<Int>(poly.size()));

        std::swap(poly[pid1], poly[pid2]);

        HashSet<Int> verts_to_update;
        verts_to_update.insert(poly[pid1].vert.begin(), poly[pid1].vert.end());
        verts_to_update.insert(poly[pid2].vert.begin(), poly[pid2].vert.end());

        HashSet<Int> edges_to_update;
        edges_to_update.insert(poly[pid1].edge.begin(), poly[pid1].edge.end());
        edges_to_update.insert(poly[pid2].edge.begin(), poly[pid2].edge.end());

        HashSet<Int> faces_to_update;
        faces_to_update.insert(poly[pid1].face.begin(), poly[pid1].face.end());
        faces_to_update.insert(poly[pid2].face.begin(), poly[pid2].face.end());

        HashSet<Int> polys_to_update;
        polys_to_update.insert(poly[pid1].poly.begin(), poly[pid1].poly.end());
        polys_to_update.insert(poly[pid2].poly.begin(), poly[pid2].poly.end());

        for (Int i : verts_to_update)
        {
            assert(0 <= i && i < static_cast<Int>(vert.size()));
            switch_data(vert[i].poly, pid1, pid2);
        }

        for (Int i : edges_to_update)
        {
            assert(0 <= i && i < static_cast<Int>(edge.size()));
            switch_data(edge[i].poly, pid1, pid2);
        }

        for (Int i : faces_to_update)
        {
            assert(0 <= i && i < static_cast<Int>(face.size()));
            switch_data(face[i].poly, pid1, pid2);
        }

        for (Int i : polys_to_update)
        {
            assert(0 <= i && i < static_cast<Int>(poly.size()));
            switch_data(poly[i].poly, pid1, pid2);
        }

        return;
    }

    template <Real TOL>
    inline Int TetMesh<TOL>::add_vert(const Point<SPACE_DIM> &v)
    {
        Int vid = find_vert(v);
        if (vid == NPOS)
        {
            vid = vert.size();
            vert.push_back(Vert{});
            vert.back().coord = v;
            vert_index[v] = vid;
        }
        return vid;
    }

    template <Real TOL>
    inline Int TetMesh<TOL>::add_edge(const Array<Int, VERTS_PER_EDGE> &vid)
    {
        assert(0 <= vid[0] && vid[0] < static_cast<Int>(vert.size()));
        assert(0 <= vid[1] && vid[1] < static_cast<Int>(vert.size()));
        assert(vid[0] != vid[1]);
        Int eid = find_edge(vid);
        if (eid == NPOS)
        {
            eid = edge.size();
            edge.push_back(Edge{});

            edge[eid].vert[0] = vid[0];
            edge[eid].vert[1] = vid[1];
            sort_by_pos<VERTS_PER_EDGE>(edge[eid].vert);

            edge[eid].edge.insert(vert[vid[0]].edge.begin(), vert[vid[0]].edge.end());
            edge[eid].edge.insert(vert[vid[1]].edge.begin(), vert[vid[1]].edge.end());
            for (Int id : edge[eid].edge)
            {
                edge[id].edge.insert(eid);
            }

            vert[vid[0]].vert.insert(vid[1]);
            vert[vid[1]].vert.insert(vid[0]);
            vert[vid[0]].edge.insert(eid);
            vert[vid[1]].edge.insert(eid);
        }
        return eid;
    }

    template <Real TOL>
    inline Int TetMesh<TOL>::add_face(const Array<Int, VERTS_PER_FACE> &vid)
    {
        assert(0 <= vid[0] && vid[0] < static_cast<Int>(vert.size()));
        assert(0 <= vid[1] && vid[1] < static_cast<Int>(vert.size()));
        assert(0 <= vid[2] && vid[2] < static_cast<Int>(vert.size()));
        assert(vid[0] != vid[1] && vid[0] != vid[2] && vid[1] != vid[2]);
        Int fid = find_face(vid);
        if (fid == NPOS)
        {
            fid = face.size();
            face.push_back(Face{});

            face[fid].vert[0] = vid[0];
            face[fid].vert[1] = vid[1];
            face[fid].vert[2] = vid[2];
            sort_by_pos<VERTS_PER_FACE>(face[fid].vert);

            const Int e01 = add_edge({vid[0], vid[1]});
            const Int e02 = add_edge({vid[0], vid[2]});
            const Int e12 = add_edge({vid[1], vid[2]});

            face[fid].edge = {e01, e02, e12};

            face[fid].face.insert(edge[e01].face.begin(), edge[e01].face.end());
            face[fid].face.insert(edge[e02].face.begin(), edge[e02].face.end());
            face[fid].face.insert(edge[e12].face.begin(), edge[e12].face.end());
            for (Int id : face[fid].face)
            {
                face[id].face.insert(fid);
            }

            vert[vid[0]].face.insert(fid);
            vert[vid[1]].face.insert(fid);
            vert[vid[2]].face.insert(fid);

            edge[e01].face.insert(fid);
            edge[e02].face.insert(fid);
            edge[e12].face.insert(fid);
        }
        return fid;
    }

    template <Real TOL>
    inline Int TetMesh<TOL>::add_poly(const Array<Int, VERTS_PER_POLY> &vid, RefineTreeNode *ptr)
    {
        assert(0 <= vid[0] && vid[0] < static_cast<Int>(vert.size()));
        assert(0 <= vid[1] && vid[1] < static_cast<Int>(vert.size()));
        assert(0 <= vid[2] && vid[2] < static_cast<Int>(vert.size()));
        assert(0 <= vid[3] && vid[3] < static_cast<Int>(vert.size()));
        assert(vid[0] != vid[1] && vid[0] != vid[2] && vid[0] != vid[3] && vid[1] != vid[2] && vid[1] != vid[3] && vid[2] != vid[3]);
        Int pid = find_poly(vid);
        if (pid == NPOS)
        {
            pid = poly.size();
            poly.push_back(Poly{});
            Int e01 = add_edge({vid[0], vid[1]});
            Int e02 = add_edge({vid[0], vid[2]});
            Int e03 = add_edge({vid[0], vid[3]});
            Int e12 = add_edge({vid[1], vid[2]});
            Int e13 = add_edge({vid[1], vid[3]});
            Int e23 = add_edge({vid[2], vid[3]});

            Int f0 = add_face({vid[1], vid[2], vid[3]});
            Int f1 = add_face({vid[0], vid[2], vid[3]});
            Int f2 = add_face({vid[0], vid[1], vid[3]});
            Int f3 = add_face({vid[0], vid[1], vid[2]});

            poly[pid].vert[0] = vid[0];
            poly[pid].vert[1] = vid[1];
            poly[pid].vert[2] = vid[2];
            poly[pid].vert[3] = vid[3];
            sort_by_pos<VERTS_PER_POLY>(poly[pid].vert);

            poly[pid].edge = {e01, e02, e03, e12, e13, e23};

            poly[pid].face = {f0, f1, f2, f3};

            poly[pid].poly.insert(face[f0].poly.begin(), face[f0].poly.end());
            poly[pid].poly.insert(face[f1].poly.begin(), face[f1].poly.end());
            poly[pid].poly.insert(face[f2].poly.begin(), face[f2].poly.end());
            poly[pid].poly.insert(face[f3].poly.begin(), face[f3].poly.end());
            for (Int id : poly[pid].poly)
            {
                poly[id].poly.insert(pid);
            }

            face[f0].poly.insert(pid);
            dangling_face.erase(f0);
            face[f1].poly.insert(pid);
            dangling_face.erase(f1);
            face[f2].poly.insert(pid);
            dangling_face.erase(f2);
            face[f3].poly.insert(pid);
            dangling_face.erase(f3);

            edge[e01].poly.insert(pid);
            dangling_edge.erase(e01);
            edge[e02].poly.insert(pid);
            dangling_edge.erase(e02);
            edge[e03].poly.insert(pid);
            dangling_edge.erase(e03);
            edge[e12].poly.insert(pid);
            dangling_edge.erase(e12);
            edge[e13].poly.insert(pid);
            dangling_edge.erase(e13);
            edge[e23].poly.insert(pid);
            dangling_edge.erase(e23);

            vert[vid[0]].poly.insert(pid);
            dangling_vert.erase(vid[0]);
            vert[vid[1]].poly.insert(pid);
            dangling_vert.erase(vid[1]);
            vert[vid[2]].poly.insert(pid);
            dangling_vert.erase(vid[2]);
            vert[vid[3]].poly.insert(pid);
            dangling_vert.erase(vid[3]);

            poly[pid].ptr = ptr;
        }
        return pid;
    }

    template <Real TOL>
    inline Int TetMesh<TOL>::add_poly(const Array<Point<SPACE_DIM>, VERTS_PER_POLY> &v, RefineTreeNode *ptr)
    {
        return add_poly({find_vert(v[0]), find_vert(v[1]), find_vert(v[2]), find_vert(v[3])}, ptr);
    }

    template <Real TOL>
    inline void TetMesh<TOL>::del_vert(const Int &vid)
    {
        assert(0 <= vid && vid < static_cast<Int>(vert.size()));
        const Int new_vid = vert.size() - 1;
        switch_vert(vid, new_vid);
        switch_data(dangling_vert, vid, new_vid);
        dangling_vert.erase(new_vid);
        const Vert &v = vert.back();
        del_edges(convert<HashSet<Int>, List<Int>>(v.edge));
        for (Int id : v.vert)
        {
            vert[id].vert.erase(new_vid);
        }
        vert_index.erase(vert[vid].coord);
        vert.pop_back();
        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::del_verts(const List<Int> &vid)
    {
        if (vid.empty())
        {
            return;
        }
        List<Int> vid_sorted = vid;
        std::ranges::sort(vid_sorted);
        vid_sorted.erase(std::ranges::unique(vid_sorted).begin(), vid_sorted.end());
        for (Int i = vid_sorted.size() - 1; i >= 0; i--)
        {
            del_vert(vid_sorted[i]);
        }
        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::del_edge(const Int &eid)
    {
        assert(0 <= eid && eid < static_cast<Int>(edge.size()));
        const Int new_eid = edge.size() - 1;
        switch_edge(eid, new_eid);
        switch_data(dangling_edge, eid, new_eid);
        dangling_edge.erase(new_eid);
        const Edge &e = edge.back();
        del_faces(convert<HashSet<Int>, List<Int>>(e.face));
        for (Int id : e.edge)
        {
            edge[id].edge.erase(new_eid);
        }
        for (Int id : e.vert)
        {
            vert[id].edge.erase(new_eid);
        }
        edge.pop_back();
        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::del_edges(const List<Int> &eid)
    {
        if (eid.empty())
        {
            return;
        }
        List<Int> eid_sorted = eid;
        std::ranges::sort(eid_sorted);
        eid_sorted.erase(std::ranges::unique(eid_sorted).begin(), eid_sorted.end());
        for (Int i = eid_sorted.size() - 1; i >= 0; i--)
        {
            del_edge(eid_sorted[i]);
        }
        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::del_face(const Int &fid)
    {
        assert(0 <= fid && fid < static_cast<Int>(face.size()));
        const Int new_fid = face.size() - 1;
        switch_face(fid, new_fid);
        switch_data(dangling_face, fid, new_fid);
        dangling_face.erase(new_fid);
        const Face &f = face.back();
        del_polys(convert<HashSet<Int>, List<Int>>(f.poly));
        for (Int id : f.face)
        {
            face[id].face.erase(new_fid);
        }
        for (Int id : f.edge)
        {
            edge[id].face.erase(new_fid);
        }
        for (Int id : f.vert)
        {
            vert[id].face.erase(new_fid);
        }
        face.pop_back();
        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::del_faces(const List<Int> &fid)
    {
        if (fid.empty())
        {
            return;
        }
        List<Int> fid_sorted = fid;
        std::ranges::sort(fid_sorted);
        fid_sorted.erase(std::ranges::unique(fid_sorted).begin(), fid_sorted.end());
        for (Int i = fid_sorted.size() - 1; i >= 0; i--)
        {
            del_face(fid_sorted[i]);
        }
        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::del_poly(const Int &pid)
    {
        assert(0 <= pid && pid < static_cast<Int>(poly.size()));
        const Int new_pid = poly.size() - 1;
        switch_poly(pid, new_pid);
        const Poly &p = poly.back();
        for (Int id : p.vert)
        {
            vert[id].poly.erase(new_pid);
            if (vert[id].poly.size() == 0)
            {
                dangling_vert.insert(id);
            }
        }
        for (Int id : p.edge)
        {
            edge[id].poly.erase(new_pid);
            if (edge[id].poly.size() == 0)
            {
                dangling_edge.insert(id);
            }
        }
        for (Int id : p.face)
        {
            face[id].poly.erase(new_pid);
            if (face[id].poly.size() == 0)
            {
                dangling_face.insert(id);
            }
        }
        for (Int id : p.poly)
        {
            poly[id].poly.erase(new_pid);
        }
        poly.pop_back();
        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::del_polys(const List<Int> &pid)
    {
        if (pid.empty())
        {
            return;
        }
        List<Int> pid_sorted = pid;
        std::ranges::sort(pid_sorted);
        pid_sorted.erase(std::ranges::unique(pid_sorted).begin(), pid_sorted.end());
        for (Int i = pid_sorted.size() - 1; i >= 0; i--)
        {
            del_poly(pid_sorted[i]);
        }
        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::clear_low_dimensional()
    {
        clear_low_dimensional_face();
        clear_low_dimensional_edge();
        clear_low_dimensional_vert();
        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::clear_low_dimensional_vert()
    {
        del_verts(convert<HashSet<Int>, List<Int>>(dangling_vert));
        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::clear_low_dimensional_edge()
    {
        del_edges(convert<HashSet<Int>, List<Int>>(dangling_edge));
        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::clear_low_dimensional_face()
    {
        del_faces(convert<HashSet<Int>, List<Int>>(dangling_face));
        return;
    }

    template <Real TOL>
    inline List<Int> TetMesh<TOL>::get_vert_opposite_edge(const Int &pid, const Int &eid)
    {
        assert(poly[pid].edge.contains(eid));
        List<Int> res;
        for (Int vid : poly[pid].vert)
        {
            if (!contain<VERTS_PER_EDGE>(edge[eid].vert, vid))
            {
                res.push_back(vid);
            }
        }
        assert(res.size() == 2);
        return res;
    }

    template <Real TOL>
    inline Int TetMesh<TOL>::get_vert_opposite_face(const Int &pid, const Int &fid)
    {
        assert(poly[pid].face.contains(fid));
        for (Int vid : poly[pid].vert)
        {
            if (!contain<VERTS_PER_FACE>(face[fid].vert, vid))
            {
                return vid;
            }
        }
        return NPOS;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::update_attributes()
    {
        for (Int i = 0; i < static_cast<Int>(vert.size()); i++)
        {
            vert[i].label = Unknown;
        }
        for (Int i = 0; i < static_cast<Int>(edge.size()); i++)
        {
            edge[i].label = Unknown;
        }
        for (Int i = 0; i < static_cast<Int>(face.size()); i++)
        {
            face[i].label = Unknown;
        }
        for (Int i = 0; i < static_cast<Int>(poly.size()); i++)
        {
            poly[i].label = Unknown;
        }
        update_face_attributes();
        update_poly_attributes();
        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::update_face_attributes()
    {
        for (Int fid = 0; fid < static_cast<Int>(face.size()); fid++)
        {
            update_face_attributes(fid);
        }
        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::update_face_attributes(const Int &fid)
    {
        List<Int> vid;
        List<Int> eid;
        List<Int> pid;
        switch (face[fid].poly.size())
        {
            case 1:
                face[fid].label = Boundary;
                vid = convert<Array<Int, VERTS_PER_FACE>, List<Int>>(face[fid].vert);
                eid = convert<HashSet<Int>, List<Int>>(face[fid].edge);
                pid = convert<HashSet<Int>, List<Int>>(face[fid].poly);
                for (Int v : vid)
                {
                    vert[v].label = Boundary;
                }
                for (Int e : eid)
                {
                    edge[e].label = Boundary;
                }
                for (Int p : pid)
                {
                    poly[p].label = Boundary;
                }
                break;
            case 2:
                face[fid].label = Interior;
                vid = convert<Array<Int, VERTS_PER_FACE>, List<Int>>(face[fid].vert);
                eid = convert<HashSet<Int>, List<Int>>(face[fid].edge);
                pid = convert<HashSet<Int>, List<Int>>(face[fid].poly);
                for (Int v : vid)
                {
                    if (vert[v].label == Unknown)
                    {
                        vert[v].label = Interior;
                    }
                }
                for (Int e : eid)
                {
                    if (edge[e].label == Unknown)
                    {
                        edge[e].label = Interior;
                    }
                }
                for (Int p : pid)
                {
                    if (poly[p].label == Unknown)
                    {
                        poly[p].label = Interior;
                    }
                }
                break;
            default:
                assert(false);
        }
        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::update_poly_attributes()
    {
        #pragma omp parallel for
        for (Int pid = 0; pid < static_cast<Int>(poly.size()); pid++)
        {
            update_poly_attributes(pid);
        }
        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::update_poly_attributes(const Int &pid)
    {
        assert(poly[pid].vert.size() == VERTS_PER_POLY);
        Array<Int, VERTS_PER_POLY> vid = poly[pid].vert;
        const Point<SPACE_DIM> &p_0 = vert[vid[0]].coord;
        const Point<SPACE_DIM> &p_1 = vert[vid[1]].coord;
        const Point<SPACE_DIM> &p_2 = vert[vid[2]].coord;
        const Point<SPACE_DIM> &p_3 = vert[vid[3]].coord;
        Matrix<Real, SPACE_DIM, SPACE_DIM> &jacobi_matrix = poly[pid].jacobi;
        jacobi_matrix(0, 0) = p_1(0) - p_0(0);
        jacobi_matrix(0, 1) = p_2(0) - p_0(0);
        jacobi_matrix(0, 2) = p_3(0) - p_0(0);
        jacobi_matrix(1, 0) = p_1(1) - p_0(1);
        jacobi_matrix(1, 1) = p_2(1) - p_0(1);
        jacobi_matrix(1, 2) = p_3(1) - p_0(1);
        jacobi_matrix(2, 0) = p_1(2) - p_0(2);
        jacobi_matrix(2, 1) = p_2(2) - p_0(2);
        jacobi_matrix(2, 2) = p_3(2) - p_0(2);
        poly[pid].inv_jacobi = jacobi_matrix.inverse();
        poly[pid].volume = std::abs(jacobi_matrix.determinant() / 6.0);
        return;
    }

    template <Real TOL>
    inline Int TetMesh<TOL>::find_vert(const Point<SPACE_DIM> &v) const
    {
        const auto it = vert_index.find(v);
        if (it != vert_index.end())
        {
            return it->second;
        }
        return NPOS;
    }

    template <Real TOL>
    inline Int TetMesh<TOL>::find_edge(const Array<Int, VERTS_PER_EDGE> &vid) const
    {
        if (vid[0] == NPOS || vid[1] == NPOS)
        {
            return NPOS;
        }
        const HashSet<Int> &list_1 = vert[vid[0]].edge;
        const HashSet<Int> &list_2 = vert[vid[1]].edge;
        for (Int eid : list_1)
        {
            if (list_2.contains(eid))
            {
                return eid;
            }
        }
        return NPOS;
    }

    template <Real TOL>
    inline Int TetMesh<TOL>::find_edge(const Array<Point<SPACE_DIM>, VERTS_PER_EDGE> &v) const
    {
        Array<Int, VERTS_PER_EDGE> vid;
        for (Int i = 0; i < VERTS_PER_EDGE; i++)
        {
            vid[i] = find_vert(v[i]);
        }
        return find_edge(vid);
    }

    template <Real TOL>
    inline Int TetMesh<TOL>::find_face(const Array<Int, VERTS_PER_FACE> &vid) const
    {
        if (vid[0] == NPOS || vid[1] == NPOS || vid[2] == NPOS)
        {
            return NPOS;
        }
        const HashSet<Int> &list_1 = vert[vid[0]].face;
        const HashSet<Int> &list_2 = vert[vid[1]].face;
        const HashSet<Int> &list_3 = vert[vid[2]].face;
        for (Int fid : list_1)
        {
            if (list_2.contains(fid) && list_3.contains(fid))
            {
                return fid;
            }
        }
        return NPOS;
    }

    template <Real TOL>
    inline Int TetMesh<TOL>::find_face(const Array<Point<SPACE_DIM>, VERTS_PER_FACE> &v) const
    {
        Array<Int, VERTS_PER_FACE> vid;
        for (Int i = 0; i < VERTS_PER_FACE; i++)
        {
            vid[i] = find_vert(v[i]);
        }
        return find_face(vid);
    }

    template <Real TOL>
    inline Int TetMesh<TOL>::find_poly(const Array<Int, VERTS_PER_POLY> &vid) const
    {
        if (vid[0] == NPOS || vid[1] == NPOS || vid[2] == NPOS || vid[3] == NPOS)
        {
            return NPOS;
        }
        const HashSet<Int> &list_1 = vert[vid[0]].poly;
        const HashSet<Int> &list_2 = vert[vid[1]].poly;
        const HashSet<Int> &list_3 = vert[vid[2]].poly;
        const HashSet<Int> &list_4 = vert[vid[3]].poly;
        for (Int pid : list_1)
        {
            if (list_2.contains(pid) && list_3.contains(pid) && list_4.contains(pid))
            {
                return pid;
            }
        }
        return NPOS;
    }

    template <Real TOL>
    inline Int TetMesh<TOL>::find_poly(const Array<Point<SPACE_DIM>, VERTS_PER_POLY> &v) const
    {
        Array<Int, VERTS_PER_POLY> vid;
        for (Int i = 0; i < VERTS_PER_POLY; i++)
        {
            vid[i] = find_vert(v[i]);
        }
        return find_poly(vid);
    }

    template <Real TOL>
    inline void TetMesh<TOL>::refine()
    {
        for (const Array<Int, VERTS_PER_POLY> &p : red_refine_poly)
        {
            Int pid = find_poly(p);
            if (pid == NPOS)
            {
                continue;
            }
            red_refine(pid);
        }
        while (!green_refine_poly.empty())
        {
            Int pid = find_poly(green_refine_poly.back());
            green_refine_poly.pop_back();
            if (pid == NPOS)
            {
                continue;
            }
            green_refine(pid);
        }
        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::red_refine(const Int &pid)
    {
        Int rpid = pid;
        if (this->poly[rpid].ptr->is_green())
        {
            rpid = recover_parent(rpid);
        }
        RefineTreeNode *p = this->poly[rpid].ptr;
        assert(p->is_red());
        for (Int tmp : neighbor(rpid))
        {
            green_refine_poly.push_back(this->poly[tmp].vert);
        }

        const List<Int> vid = convert<Array<Int, VERTS_PER_POLY>, List<Int>>(this->poly[rpid].vert);
        const Int v_0 = vid[0];
        const Int v_1 = vid[1];
        const Int v_2 = vid[2];
        const Int v_3 = vid[3];
        const Int v_01 = this->add_vert((this->vert[v_0].coord + this->vert[v_1].coord) / 2.0);
        const Int v_02 = this->add_vert((this->vert[v_0].coord + this->vert[v_2].coord) / 2.0);
        const Int v_03 = this->add_vert((this->vert[v_0].coord + this->vert[v_3].coord) / 2.0);
        const Int v_12 = this->add_vert((this->vert[v_1].coord + this->vert[v_2].coord) / 2.0);
        const Int v_13 = this->add_vert((this->vert[v_1].coord + this->vert[v_3].coord) / 2.0);
        const Int v_23 = this->add_vert((this->vert[v_2].coord + this->vert[v_3].coord) / 2.0);

        const Array<Int, 6> ov = {v_01, v_02, v_03, v_12, v_13, v_23};

        const Array<Real, 3> len = {
                    (this->vert[v_01].coord - this->vert[v_23].coord).template lpNorm<2>(),
                    (this->vert[v_02].coord - this->vert[v_13].coord).template lpNorm<2>(),
                    (this->vert[v_03].coord - this->vert[v_12].coord).template lpNorm<2>()
                };

        constexpr Int pair[3][2] = {{0, 5}, {1, 4}, {2, 3}};

        Int l = 0, s = 0, m = 0;
        for (Int i = 0; i < 3; i++)
        {
            if (len[l] < len[i])
            {
                l = i;
            }
            if (len[s] >= len[i])
            {
                s = i;
            }
        }
        m = 3 - l - s;
        assert(l != s);
        Array<Int, 6> v = {ov[pair[s][0]], ov[pair[m][0]], ov[pair[s][1]], ov[pair[m][1]], ov[pair[l][0]], ov[pair[l][1]]};

        p->set_child({
                    {this->vert[v_0].coord, this->vert[v_01].coord, this->vert[v_02].coord, this->vert[v_03].coord},
                    {this->vert[v_1].coord, this->vert[v_01].coord, this->vert[v_12].coord, this->vert[v_13].coord},
                    {this->vert[v_2].coord, this->vert[v_02].coord, this->vert[v_12].coord, this->vert[v_23].coord},
                    {this->vert[v_3].coord, this->vert[v_03].coord, this->vert[v_13].coord, this->vert[v_23].coord},

                    {this->vert[v[0]].coord, this->vert[v[1]].coord, this->vert[v[2]].coord, this->vert[v[4]].coord},
                    {this->vert[v[0]].coord, this->vert[v[1]].coord, this->vert[v[2]].coord, this->vert[v[5]].coord},
                    {this->vert[v[0]].coord, this->vert[v[2]].coord, this->vert[v[3]].coord, this->vert[v[4]].coord},
                    {this->vert[v[0]].coord, this->vert[v[2]].coord, this->vert[v[3]].coord, this->vert[v[5]].coord}
                });

        this->add_poly({v_0, v_01, v_02, v_03}, p->child[0]);
        this->add_poly({v_1, v_01, v_12, v_13}, p->child[1]);
        this->add_poly({v_2, v_02, v_12, v_23}, p->child[2]);
        this->add_poly({v_3, v_03, v_13, v_23}, p->child[3]);
        this->add_poly({v[0], v[1], v[2], v[4]}, p->child[4]);
        this->add_poly({v[0], v[1], v[2], v[5]}, p->child[5]);
        this->add_poly({v[0], v[2], v[3], v[4]}, p->child[6]);
        this->add_poly({v[0], v[2], v[3], v[5]}, p->child[7]);

        this->del_poly(rpid);

        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::green_refine(const Int &pid)
    {
        Int rpid = pid;
        if (this->poly[rpid].ptr->is_green())
        {
            rpid = recover_parent(rpid);
        }

        constexpr Int V_01 = 1;
        constexpr Int V_02 = 10;
        constexpr Int V_03 = 100;
        constexpr Int V_12 = 1000;
        constexpr Int V_13 = 10000;
        constexpr Int V_23 = 100000;

        List<Int> vid = convert<Array<Int, VERTS_PER_POLY>, List<Int>>(this->poly[rpid].vert);
        List<Point<SPACE_DIM>> v;
        for (Int i : vid)
        {
            v.push_back(this->vert[i].coord);
        }

        v.push_back((v[0] + v[1]) / 2.0);
        vid.push_back(this->find_vert(v.back()));
        const Int v_01 = vid.back() == NPOS ? 0 : V_01;
        v.push_back((v[0] + v[2]) / 2.0);
        vid.push_back(this->find_vert(v.back()));
        const Int v_02 = vid.back() == NPOS ? 0 : V_02;
        v.push_back((v[0] + v[3]) / 2.0);
        vid.push_back(this->find_vert(v.back()));
        const Int v_03 = vid.back() == NPOS ? 0 : V_03;
        v.push_back((v[1] + v[2]) / 2.0);
        vid.push_back(this->find_vert(v.back()));
        const Int v_12 = vid.back() == NPOS ? 0 : V_12;
        v.push_back((v[1] + v[3]) / 2.0);
        vid.push_back(this->find_vert(v.back()));
        const Int v_13 = vid.back() == NPOS ? 0 : V_13;
        v.push_back((v[2] + v[3]) / 2.0);
        vid.push_back(this->find_vert(v.back()));
        const Int v_23 = vid.back() == NPOS ? 0 : V_23;

        if ((v_01 + v_02 + v_03 + v_12 + v_13 + v_23) == 0)
        {
            return;
        }

        List<Int> list;
        for (Int k = 0; k < static_cast<Int>(TET_TWICE_REFINE_TABLE.size()); k++)
        {
            Point<SPACE_DIM> p = Point<SPACE_DIM>::Zero();
            for (Int i = 0; i < 4; i++)
            {
                p += (TET_TWICE_REFINE_TABLE[k].first[i] * this->vert[vid[i]].coord);
            }
            if (this->find_vert(p) != NPOS)
            {
                list.push_back(k);
            }
        }
        if (!list.empty())
        {
            red_refine(rpid);
            for (Int i = 0; i < static_cast<Int>(vid.size()); i++)
            {
                if (vid[i] == NPOS)
                {
                    vid[i] = this->find_vert(v[i]);
                }
            }
            for (Int i : list)
            {
                Int eid = this->find_edge({vid[TET_TWICE_REFINE_TABLE[i].second[0]], vid[TET_TWICE_REFINE_TABLE[i].second[1]]});
                for (Int npid : this->edge[eid].poly)
                {
                    green_refine_poly.push_back(this->poly[npid].vert);
                }
            }
            return;
        }

        switch (v_01 + v_02 + v_03 + v_12 + v_13 + v_23)
        {
            case V_01:
                green_refine_edge(rpid, this->find_edge({vid[0], vid[1]}));
                break;
            case V_02:
                green_refine_edge(rpid, this->find_edge({vid[0], vid[2]}));
                break;
            case V_03:
                green_refine_edge(rpid, this->find_edge({vid[0], vid[3]}));
                break;
            case V_12:
                green_refine_edge(rpid, this->find_edge({vid[1], vid[2]}));
                break;
            case V_13:
                green_refine_edge(rpid, this->find_edge({vid[1], vid[3]}));
                break;
            case V_23:
                green_refine_edge(rpid, this->find_edge({vid[2], vid[3]}));
                break;

            case V_01 + V_02:
                green_refine_face(rpid, this->find_face({vid[0], vid[1], vid[2]}));
                break;
            case V_01 + V_03:
                green_refine_face(rpid, this->find_face({vid[0], vid[1], vid[3]}));
                break;
            case V_02 + V_03:
                green_refine_face(rpid, this->find_face({vid[0], vid[2], vid[3]}));
                break;
            case V_12 + V_13:
                green_refine_face(rpid, this->find_face({vid[1], vid[2], vid[3]}));
                break;

            case V_01 + V_12:
                green_refine_face(rpid, this->find_face({vid[0], vid[1], vid[2]}));
                break;
            case V_01 + V_13:
                green_refine_face(rpid, this->find_face({vid[0], vid[1], vid[3]}));
                break;
            case V_02 + V_23:
                green_refine_face(rpid, this->find_face({vid[0], vid[2], vid[3]}));
                break;
            case V_12 + V_23:
                green_refine_face(rpid, this->find_face({vid[1], vid[2], vid[3]}));
                break;

            case V_02 + V_12:
                green_refine_face(rpid, this->find_face({vid[0], vid[1], vid[2]}));
                break;
            case V_03 + V_13:
                green_refine_face(rpid, this->find_face({vid[0], vid[1], vid[3]}));
                break;
            case V_03 + V_23:
                green_refine_face(rpid, this->find_face({vid[0], vid[2], vid[3]}));
                break;
            case V_13 + V_23:
                green_refine_face(rpid, this->find_face({vid[1], vid[2], vid[3]}));
                break;

            case V_01 + V_02 + V_12:
                green_refine_face(rpid, this->find_face({vid[0], vid[1], vid[2]}));
                break;
            case V_01 + V_03 + V_13:
                green_refine_face(rpid, this->find_face({vid[0], vid[1], vid[3]}));
                break;
            case V_02 + V_03 + V_23:
                green_refine_face(rpid, this->find_face({vid[0], vid[2], vid[3]}));
                break;
            case V_12 + V_13 + V_23:
                green_refine_face(rpid, this->find_face({vid[1], vid[2], vid[3]}));
                break;

            default:
                red_refine(rpid);
                break;
        }

        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::green_refine_edge(const Int &pid, const Int &eid)
    {
        RefineTreeNode *p = this->poly[pid].ptr;
        assert(p->is_red());
        List<Int> vid = convert<Array<Int, VERTS_PER_EDGE>, List<Int>>(this->edge[eid].vert);
        const Int v_0 = vid[0];
        const Int v_1 = vid[1];
        vid = this->get_vert_opposite_edge(pid, eid);
        const Int v_2 = vid[0];
        const Int v_3 = vid[1];
        const Int v_01 = this->add_vert((this->vert[v_0].coord + this->vert[v_1].coord) / 2.0);
        p->set_child({
                    {this->vert[v_2].coord, this->vert[v_3].coord, this->vert[v_0].coord, this->vert[v_01].coord},
                    {this->vert[v_2].coord, this->vert[v_3].coord, this->vert[v_1].coord, this->vert[v_01].coord}
                });
        this->add_poly({v_2, v_3, v_0, v_01}, p->child[0]);
        this->add_poly({v_2, v_3, v_1, v_01}, p->child[1]);
        this->del_poly(pid);
        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::green_refine_face(const Int &pid, const Int &fid)
    {
        RefineTreeNode *p = this->poly[pid].ptr;
        assert(p->is_red());
        List<Int> vid = convert<Array<Int, VERTS_PER_FACE>, List<Int>>(this->face[fid].vert);
        const Int v_0 = vid[0];
        const Int v_1 = vid[1];
        const Int v_2 = vid[2];
        const Int v_3 = this->get_vert_opposite_face(pid, fid);
        const Int v_01 = this->add_vert((this->vert[v_0].coord + this->vert[v_1].coord) / 2.0);
        const Int v_02 = this->add_vert((this->vert[v_0].coord + this->vert[v_2].coord) / 2.0);
        const Int v_12 = this->add_vert((this->vert[v_1].coord + this->vert[v_2].coord) / 2.0);
        p->set_child({
                    {this->vert[v_3].coord, this->vert[v_0].coord, this->vert[v_01].coord, this->vert[v_02].coord},
                    {this->vert[v_3].coord, this->vert[v_1].coord, this->vert[v_01].coord, this->vert[v_12].coord},
                    {this->vert[v_3].coord, this->vert[v_2].coord, this->vert[v_02].coord, this->vert[v_12].coord},
                    {this->vert[v_3].coord, this->vert[v_01].coord, this->vert[v_02].coord, this->vert[v_12].coord}
                });
        this->add_poly({v_3, v_0, v_01, v_02}, p->child[0]);
        this->add_poly({v_3, v_1, v_01, v_12}, p->child[1]);
        this->add_poly({v_3, v_2, v_02, v_12}, p->child[2]);
        this->add_poly({v_3, v_01, v_02, v_12}, p->child[3]);
        this->del_poly(pid);
        for (Int eid : this->face[fid].edge)
        {
            for (Int npid : this->edge[eid].poly)
            {
                green_refine_poly.push_back(this->poly[npid].vert);
            }
        }
        return;
    }

    template <Real TOL>
    inline Int TetMesh<TOL>::recover_parent(const Int &pid)
    {
        assert(this->poly[pid].ptr->parent != nullptr);
        RefineTreeNode *parent = this->poly[pid].ptr->parent;
        for (const RefineTreeNode *p : parent->child)
        {
            this->del_poly(this->find_poly(p->data));
        }
        parent->del_child();
        return add_poly(parent->data, parent);
    }

    template <Real TOL>
    inline List<Int> TetMesh<TOL>::neighbor(const Int &pid)
    {
        HashSet<Int> res = this->poly[pid].poly;
        for (Int eid : this->poly[pid].edge)
        {
            res.insert(this->edge[eid].poly.begin(), this->edge[eid].poly.end());
        }
        res.erase(pid);
        return convert<HashSet<Int>, List<Int>>(res);
    }

    template <Real TOL>
    inline void TetMesh<TOL>::build_tree_mapping()
    {
        List<RefineTreeNode *> list = tree.root;
        while (!list.empty())
        {
            RefineTreeNode *ptr = list.back();
            list.pop_back();
            const Int pid = this->find_poly(ptr->data);
            if (pid != NPOS)
            {
                this->poly[pid].ptr = ptr;
            }
            list.insert(list.end(), ptr->child.begin(), ptr->child.end());
        }
    }

    template <Real TOL>
    inline TetMesh<TOL>::TetMesh(const char *filename)
    {
        load(filename);
    }

    template <Real TOL>
    inline TetMesh<TOL>::TetMesh(const String &filename)
    {
        load(filename);
    }

    template <Real TOL>
    inline TetMesh<TOL>::TetMesh(const TetMesh &mesh)
    {
        dangling_face = mesh.dangling_face;
        dangling_edge = mesh.dangling_edge;
        dangling_vert = mesh.dangling_vert;
        vert_index = mesh.vert_index;
        vert = mesh.vert;
        edge = mesh.edge;
        face = mesh.face;
        poly = mesh.poly;
        tree = mesh.tree;
        red_refine_poly = mesh.red_refine_poly;
        green_refine_poly = mesh.green_refine_poly;
        build_tree_mapping();
    }

    template <Real TOL>
    inline void TetMesh<TOL>::load(const String &filename)
    {
        if (const String ext = to_lower(get_file_extension(filename)); ext == "mesh")
        {
            load_mesh(filename);
        }
        else
        {
            ASN_ERROR("void TetMesh<TOL>::load() : File format not supported yet.");
        }
        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::save(const String &filename) const
    {
        if (const String ext = to_lower(get_file_extension(filename)); ext == "mesh")
        {
            save_mesh(filename);
        }
        else if (ext == "obj")
        {
            save_obj(filename);
        }
        else
        {
            ASN_ERROR("void TetMesh<TOL>::save(const String &filename) const : File format not supported yet.");
        }
        return;
    }

    template <Real TOL>
    inline Real TetMesh<TOL>::get_face_diameter(const Int &pid) const
    {
        assert(0 <= pid && pid < static_cast<Int>(this->face.size()));
        List<const Point<SPACE_DIM> *> point;
        for (const Int vid : this->face[pid].vert)
        {
            point.push_back(&this->vert[vid].coord);
        }
        Real diam = 0.0;
        for (Int i = 0; i < static_cast<Int>(point.size()); i++)
        {
            for (Int j = i + 1; j < static_cast<Int>(point.size()); j++)
            {
                const Real tmp = (*point[i] - *point[j]).lpNorm<2>();
                if (tmp >= diam)
                {
                    diam = tmp;
                }
            }
        }
        return diam;
    }

    template <Real TOL>
    inline Real TetMesh<TOL>::get_poly_diameter(const Int &pid) const
    {
        assert(0 <= pid && pid < static_cast<Int>(this->poly.size()));
        List<const Point<SPACE_DIM> *> point;
        for (const Int vid : this->poly[pid].vert)
        {
            point.push_back(&this->vert[vid].coord);
        }
        Real diam = 0.0;
        for (Int i = 0; i < static_cast<Int>(point.size()); i++)
        {
            for (Int j = i + 1; j < static_cast<Int>(point.size()); j++)
            {
                const Real tmp = (*point[i] - *point[j]).lpNorm<2>();
                if (tmp >= diam)
                {
                    diam = tmp;
                }
            }
        }
        return diam;
    }

    template <Real TOL>
    inline Int TetMesh<TOL>::num_vert() const
    {
        return vert.size();
    }

    template <Real TOL>
    inline Int TetMesh<TOL>::num_edge() const
    {
        return edge.size();
    }

    template <Real TOL>
    inline Int TetMesh<TOL>::num_face() const
    {
        return face.size();
    }

    template <Real TOL>
    inline Int TetMesh<TOL>::num_poly() const
    {
        return poly.size();
    }

    template <Real TOL>
    inline const TetMesh<TOL>::Vert &TetMesh<TOL>::get_vert(const Int &vid) const
    {
        assert(0 <= vid && vid < static_cast<Int>(vert.size()));
        return vert[vid];
    }

    template <Real TOL>
    inline const TetMesh<TOL>::Edge &TetMesh<TOL>::get_edge(const Int &eid) const
    {
        assert(0 <= eid && eid < static_cast<Int>(edge.size()));
        return edge[eid];
    }

    template <Real TOL>
    inline const TetMesh<TOL>::Face &TetMesh<TOL>::get_face(const Int &fid) const
    {
        assert(0 <= fid && fid < static_cast<Int>(face.size()));
        return face[fid];
    }

    template <Real TOL>
    inline const TetMesh<TOL>::Poly &TetMesh<TOL>::get_poly(const Int &pid) const
    {
        assert(0 <= pid && pid < static_cast<Int>(poly.size()));
        return poly[pid];
    }

    template <Real TOL>
    inline Bool TetMesh<TOL>::check_regular(const Int &component, const Int &cavity, const Int &hole)
    {
        std::cout << "Checking Mesh Topology ..." << std::endl;

        update_face_attributes();

        Int mesh_component = 0;
        List<Bool> label(vert.size(), true);
        List<Int> queue;
        for (Int i = 0; i < static_cast<Int>(label.size()); i++)
        {
            if (label[i])
            {
                mesh_component++;
                queue.push_back(i);
                while (!queue.empty())
                {
                    Int vid = queue.back();
                    queue.pop_back();

                    label[vid] = false;

                    for (Int v : vert[vid].vert)
                    {
                        if (label[v])
                        {
                            queue.push_back(v);
                        }
                    }
                }
            }
        }

        std::cout << std::left << std::setw(40) << "    Num of mesh component = " << mesh_component << " / " << component << "." << std::endl;

        Int boundary_component = 0;
        for (Int i = 0; i < static_cast<Int>(label.size()); i++)
        {
            if (vert[i].label == Boundary)
            {
                label[i] = true;
            }
        }
        queue.clear();
        for (Int i = 0; i < static_cast<Int>(label.size()); i++)
        {
            if (label[i])
            {
                boundary_component++;
                queue.push_back(i);
                while (!queue.empty())
                {
                    Int vid = queue.back();
                    queue.pop_back();

                    label[vid] = false;

                    for (Int v : vert[vid].vert)
                    {
                        if (label[v])
                        {
                            queue.push_back(v);
                        }
                    }
                }
            }
        }
        std::cout << std::left << std::setw(40) << "    Num of boundary component = " << boundary_component << " / " << (component + cavity) << "." << std::endl;

        Int euler_characteristic = vert.size() - edge.size() + face.size() - poly.size();

        std::cout << std::left << std::setw(40) << "    Num of cavity = " << (boundary_component - mesh_component) << " / " << cavity << "." << std::endl;
        std::cout << std::left << std::setw(40) << "    Num of hole = " << (boundary_component - euler_characteristic) << " / " << hole << "." << std::endl;
        std::cout << std::left << std::setw(40) << "    Euler characteristic = " << euler_characteristic << " / " << (component - hole + cavity) << "." << std::endl;

        Bool is_regular = true;
        if (component != mesh_component)
        {
            is_regular = false;
        }
        if (euler_characteristic != (component + cavity - hole))
        {
            is_regular = false;
        }
        return is_regular;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::init_h_adaptive()
    {
        for (Int pid = 0; pid < static_cast<Int>(this->poly.size()); pid++)
        {
            Array<Int, VERTS_PER_POLY> vid = this->poly[pid].vert;
            this->poly[pid].ptr = tree.add_root_node({this->vert[vid[0]].coord, this->vert[vid[1]].coord, this->vert[vid[2]].coord, this->vert[vid[3]].coord});
        }
        clear_h_adaptive_flag();
        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::clear_h_adaptive_flag()
    {
        #pragma omp parallel for
        for (Int pid = 0; pid < static_cast<Int>(this->poly.size()); pid++)
        {
            this->poly[pid].refine = false;
            this->poly[pid].coarsen = false;
        }
        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::set_refine_flag(const Int &pid)
    {
        this->poly[pid].refine = true;
        return;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::set_coarsen_flag(const Int &pid)
    {
        this->poly[pid].coarsen = true;
        return;
    }

    template <Real TOL>
    inline Bool TetMesh<TOL>::h_adaptive()
    {
        Bool changed = false;
        red_refine_poly.clear();
        green_refine_poly.clear();
        for (Int pid = 0; pid < static_cast<Int>(this->poly.size()); pid++)
        {
            if (this->poly[pid].refine)
            {
                changed = true;
                red_refine_poly.push_back(this->poly[pid].vert);
            }
        }
        refine();
        clear_h_adaptive_flag();
        this->clear_low_dimensional();
        std::cout << "Current mesh: [ " << this->vert.size() << " V / ";
        std::cout << this->edge.size() << " E / ";
        std::cout << this->face.size() << " F / ";
        std::cout << this->poly.size() << " P ] " << std::endl;
        return changed;
    }

    template <Real TOL>
    inline void TetMesh<TOL>::global_refine()
    {
        red_refine_poly.clear();
        green_refine_poly.clear();
        for (Int pid = 0; pid < static_cast<Int>(this->poly.size()); pid++)
        {
            red_refine_poly.push_back(this->poly[pid].vert);
        }
        refine();
        clear_h_adaptive_flag();
        this->clear_low_dimensional();
        std::cout << "Current mesh: [ " << this->vert.size() << " V / ";
        std::cout << this->edge.size() << " E / ";
        std::cout << this->face.size() << " F / ";
        std::cout << this->poly.size() << " P ] " << std::endl;
        return;
    }

    template <Real TOL>
    inline TetMesh<TOL> &TetMesh<TOL>::operator=(const TetMesh &mesh)
    {
        if (this != &mesh)
        {
            dangling_face = mesh.dangling_face;
            dangling_edge = mesh.dangling_edge;
            dangling_vert = mesh.dangling_vert;
            vert_index = mesh.vert_index;
            vert = mesh.vert;
            edge = mesh.edge;
            face = mesh.face;
            poly = mesh.poly;
            tree = mesh.tree;
            red_refine_poly = mesh.red_refine_poly;
            green_refine_poly = mesh.green_refine_poly;
            build_tree_mapping();
        }
        return *this;
    }
};

#endif
