#ifndef ASN_MESH_TETMESH_H
#define ASN_MESH_TETMESH_H

#include "Utils.h"
#include "Tree.h"

namespace Asn::Mesh
{
    template <Real TOL = DEFAULT_REL_TOL>
    class TetMesh;

    template <Real TOL>
    class TetMesh
    {
    public:
        struct Vert;
        struct Edge;
        struct Face;
        struct Poly;

    public:
        static constexpr Int MESH_DIM = 3;
        static constexpr Int SPACE_DIM = 3;

        static constexpr Int VERTS_PER_EDGE = 2;
        static constexpr Int VERTS_PER_FACE = 3;
        static constexpr Int EDGES_PER_FACE = 3;
        static constexpr Int VERTS_PER_POLY = 4;
        static constexpr Int EDGES_PER_POLY = 6;
        static constexpr Int FACES_PER_POLY = 4;

        static constexpr Int NUM_RED_CHILD = 8;

    protected:
        static const CmpPoint<SPACE_DIM, TOL> CMP_POINT;

    public:
        using RefineTree = Tree<Array<Point<SPACE_DIM>, 4>, NUM_RED_CHILD>;
        using RefineTreeNode = Tree<Array<Point<SPACE_DIM>, 4>, NUM_RED_CHILD>::Node;

    public:
        struct Vert
        {
            Point<SPACE_DIM> coord = Point<SPACE_DIM>::Zero();

            HashSet<Int> vert;
            HashSet<Int> edge;
            HashSet<Int> face;
            HashSet<Int> poly;

            Int label = Unknown;
        };

        struct Edge
        {
            Array<Int, VERTS_PER_EDGE> vert;
            HashSet<Int> edge;
            HashSet<Int> face;
            HashSet<Int> poly;

            Real length;

            Int label = Unknown;
        };

        struct Face
        {
            Array<Int, VERTS_PER_FACE> vert;
            HashSet<Int> edge;
            HashSet<Int> face;
            HashSet<Int> poly;

            Real area;

            Int label = Unknown;
        };

        struct Poly
        {
            Array<Int, VERTS_PER_POLY> vert;
            HashSet<Int> edge;
            HashSet<Int> face;
            HashSet<Int> poly;

            Int label = Unknown;

            Real volume;
            Matrix<Real, SPACE_DIM, SPACE_DIM> jacobi = Matrix<Real, SPACE_DIM, SPACE_DIM>::Zero();
            Matrix<Real, SPACE_DIM, SPACE_DIM> inv_jacobi = Matrix<Real, SPACE_DIM, SPACE_DIM>::Zero();

            Bool refine = false;
            Bool coarsen = false;
            RefineTreeNode *ptr = nullptr;
        };

    protected:
        HashSet<Int> dangling_face;
        HashSet<Int> dangling_edge;
        HashSet<Int> dangling_vert;

        Map<Point<SPACE_DIM>, Int, CmpPoint<SPACE_DIM, TOL>> vert_index;

        List<Vert> vert;
        List<Edge> edge;
        List<Face> face;
        List<Poly> poly;

        RefineTree tree;

        List<Array<Int, VERTS_PER_POLY>> red_refine_poly;
        List<Array<Int, VERTS_PER_POLY>> green_refine_poly;

    protected:
        Order cmp_vert_idx(const Int &i, const Int &j) const;
        template <Int N>
        void sort_by_pos(Array<Int, N> &list) const;

        void load_mesh(const String &filename);
        void save_mesh(const String &filename) const;
        void save_obj(const String &filename) const;

        void clear();

        void switch_vert(const Int &vid1, const Int &vid2);
        void switch_edge(const Int &eid1, const Int &eid2);
        void switch_face(const Int &fid1, const Int &fid2);
        void switch_poly(const Int &pid1, const Int &pid2);

        Int add_vert(const Point<SPACE_DIM> &v);
        Int add_edge(const Array<Int, VERTS_PER_EDGE> &vid);
        Int add_face(const Array<Int, VERTS_PER_FACE> &vid);
        Int add_poly(const Array<Int, VERTS_PER_POLY> &vid, RefineTreeNode *ptr);
        Int add_poly(const Array<Point<SPACE_DIM>, VERTS_PER_POLY> &v, RefineTreeNode *ptr);

        void del_vert(const Int &vid);
        void del_verts(const List<Int> &vid);
        void del_edge(const Int &eid);
        void del_edges(const List<Int> &eid);
        void del_face(const Int &fid);
        void del_faces(const List<Int> &fid);
        void del_poly(const Int &pid);
        void del_polys(const List<Int> &pid);

        void clear_low_dimensional();
        void clear_low_dimensional_vert();
        void clear_low_dimensional_edge();
        void clear_low_dimensional_face();

        List<Int> get_vert_opposite_edge(const Int &pid, const Int &eid);
        Int get_vert_opposite_face(const Int &pid, const Int &fid);

        void update_attributes();
        void update_face_attributes();
        void update_face_attributes(const Int &fid);
        void update_poly_attributes();
        void update_poly_attributes(const Int &pid);

        Int find_vert(const Point<SPACE_DIM> &v) const;
        Int find_edge(const Array<Int, VERTS_PER_EDGE> &vid) const;
        Int find_edge(const Array<Point<SPACE_DIM>, VERTS_PER_EDGE> &v) const;
        Int find_face(const Array<Int, VERTS_PER_FACE> &vid) const;
        Int find_face(const Array<Point<SPACE_DIM>, VERTS_PER_FACE> &v) const;
        Int find_poly(const Array<Int, VERTS_PER_POLY> &vid) const;
        Int find_poly(const Array<Point<SPACE_DIM>, VERTS_PER_POLY> &v) const;

        void refine();
        void red_refine(const Int &pid);
        void green_refine(const Int &pid);
        void green_refine_edge(const Int &pid, const Int &eid);
        void green_refine_face(const Int &pid, const Int &fid);
        Int recover_parent(const Int &pid);
        List<Int> neighbor(const Int &pid);

        void build_tree_mapping();

    public:
        TetMesh() = default;
        explicit TetMesh(const char *filename);
        explicit TetMesh(const String &filename);
        TetMesh(const TetMesh &mesh);
        virtual ~TetMesh() = default;

        void load(const String &filename);
        void save(const String &filename) const;

        Real get_face_diameter(const Int &pid) const;
        Real get_poly_diameter(const Int &pid) const;

        Int num_vert() const;
        Int num_edge() const;
        Int num_face() const;
        Int num_poly() const;

        const Vert &get_vert(const Int &vid) const;
        const Edge &get_edge(const Int &eid) const;
        const Face &get_face(const Int &fid) const;
        const Poly &get_poly(const Int &pid) const;

        Bool check_regular(const Int &component = 1, const Int &cavity = 0, const Int &hole = 0);

        void init_h_adaptive();
        void clear_h_adaptive_flag();
        void set_refine_flag(const Int &pid);
        void set_coarsen_flag(const Int &pid);
        Bool h_adaptive();
        void global_refine();

        TetMesh &operator=(const TetMesh &mesh);
    };
};

#endif
