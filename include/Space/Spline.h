#ifndef ASN_MATH_SEM_SPLINE_H
#define ASN_MATH_SEM_SPLINE_H

#include "Space.h"

namespace Asn::Math
{
    using BaseMesh = Mesh::TetMesh<>;
    using ConstNodePtr = const BaseMesh::RefineTreeNode *;

    template <Real TOL = DEFAULT_REL_TOL>
    class Spline;

    template <Real TOL>
    class Spline : public BaseMesh
    {
    public:
        static Bool inside_poly(const Point<BaseMesh::MESH_DIM> &p, const Array<Point<BaseMesh::MESH_DIM>, 4> &t);
        static Bool inside_poly(const List<Point<BaseMesh::MESH_DIM>> &p, const Array<Point<BaseMesh::MESH_DIM>, 4> &t);

    private:
        Int max_order;

        List<Index<2>> index_face;
        List<Index<3>> index_poly;
        List<Index<3>> index_poly_all;

        List<Int> num_dof_vert;
        List<List<Int>> num_dof_edge;
        List<IndexedData<2, Int>> num_dof_face;
        List<IndexedData<3, Int>> num_dof_poly;
        List<List<Int>> num_dof_poly_all;

        Matrix<Real> solution;

    public:
        Spline() = default;
        Spline(const Spline &spline) = default;
        virtual ~Spline() override = default;

        void init(const Space<TOL> &space);
        void init(const Matrix<Real> &solution);

        Real get_value(const Point<BaseMesh::MESH_DIM> &p, const Int &index = 0, ConstNodePtr ptr = nullptr) const;
        Vector<Real> get_value(const Point<BaseMesh::MESH_DIM> &p, ConstNodePtr ptr = nullptr) const;

        ConstNodePtr find_node(const List<Point<BaseMesh::MESH_DIM>> &p, ConstNodePtr ptr = nullptr) const;
        ConstNodePtr find_leaf(const List<Point<BaseMesh::MESH_DIM>> &p, ConstNodePtr ptr = nullptr) const;

        Matrix<Real> interp(const Matrix<Real> &solution, const Space<TOL> &space);

        Spline &operator=(const Spline &spline) = delete;
    };
};

#include "Spline.hpp"

#endif
