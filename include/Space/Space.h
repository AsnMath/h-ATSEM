#ifndef ASN_MATH_SPACE_SPACE_H
#define ASN_MATH_SPACE_SPACE_H

#include "Index.h"

#include "../Config.h"
#include "../Quadrature/Quadrature.h"
#include "../Mesh/TetMesh.h"

namespace Asn::Math
{
    using BaseMesh = Mesh::TetMesh<>;

    template <Real TOL = DEFAULT_REL_TOL>
    class Space;

    template <Real TOL>
    class Space final : public Mesh::TetMesh<TOL>
    {
    public:
        template <typename Type>
        static Type jacobi(const Int &a, const Int &b, const Int &k, const Type &tu, const Type &td);
        template <typename Type>
        static Type koornwinder(const Index<3> &i, const Point<3, Type> &x);
        template <typename Type>
        static Type psi(const Index<3> &i, const Point<3, Type> &x);
        static Point<3> grad_psi(const Index<3> &i, const Point<3> &x);
        static Pair<Point<3>, Matrix<Real, 3, 3>> grad_hess_psi(const Index<3> &i, const Point<3> &x);

        static List<Index<2>> get_index_list_face(const Int &order);
        static List<Index<3>> get_index_list_poly(const Int &order);
        static List<Index<3>> get_index_list_poly_all(const Int &order);

        static SparseMatrix<Real> map_matrix(const SparseMatrix<Real> &A, const HashMap<Int, Int> &index, const Int &size = 0);
        static Matrix<Real> map_matrix_row(const Matrix<Real> &A, const HashMap<Int, Int> &index, const Int &row = 0);

    public:
        const Int MAX_ORDER;
        const Int QUAD_ORDER;
        const Quadrature QUAD;

        const List<Index<2>> INDEX_FACE;
        const List<Index<3>> INDEX_POLY;
        const List<Index<3>> INDEX_POLY_ALL;

    protected:
        Int _num_dof;

        List<Int> _num_dof_vert;
        List<List<Int>> _num_dof_edge;
        List<IndexedData<2, Int>> _num_dof_face;
        List<IndexedData<3, Int>> _num_dof_poly;
        List<List<Int>> _num_dof_poly_all;

        List<List<Real>> _inner_product_value;
        IndexedData<3, List<Real>> _psi_value;
        IndexedData<3, List<Point<3>>> _grad_psi_value;
        IndexedData<3, Array<List<Point<3>>, 4>> _grad_psi_value_face;
        IndexedData<3, List<Matrix<Real, 3, 3>>> _hess_psi_value;

        List<IndexedData<3, List<Point<3>>>> _trans_grad_psi_value;

        List<Int> _boundary_vert;
        List<Int> _boundary_edge;
        List<Int> _boundary_face;
        HashSet<Int> _boundary_dof;

        SparseMatrix<Real> _stiff_matrix;
        SparseMatrix<Real> _mass_matrix;

        List<List<Real>> jacobi_edge;
        List<List<Real>> jacobi_face_1;
        List<List<Real>> jacobi_face_2;

    public:
        const Int &num_dof;

        const List<Int> &num_dof_vert;
        const List<List<Int>> &num_dof_edge;
        const List<IndexedData<2, Int>> &num_dof_face;
        const List<IndexedData<3, Int>> &num_dof_poly;
        const List<List<Int>> &num_dof_poly_all;

        const List<List<Real>> &inner_product_value;
        const IndexedData<3, List<Real>> &psi_value;
        const IndexedData<3, List<Point<3>>> &grad_psi_value;
        const IndexedData<3, Array<List<Point<3>>, 4>> &grad_psi_value_face;
        const IndexedData<3, List<Matrix<Real, 3, 3>>> &hess_psi_value;

        const List<IndexedData<3, List<Point<3>>>> &trans_grad_psi_value;

        const List<Int> &boundary_vert;
        const List<Int> &boundary_edge;
        const List<Int> &boundary_face;
        const HashSet<Int> &boundary_dof;

        const SparseMatrix<Real> &stiff_matrix;
        const SparseMatrix<Real> &mass_matrix;

    protected:
        void build_dof();
        void build_precomputation();
        void build_precomputation_element();
        void build_boundary();

    public:
        Space(const String &mesh_file, const Int &order, const Int &quad_order);
        Space(const Space &space);
        virtual ~Space() override = default;

        void init();

        Real build_element_stiff_matrix(const Int &pid, const Index<3> &i, const Index<3> &j) const;
        void build_stiff_matrix();
        void build_mass_matrix();

        template <typename VectorLike>
        void interp_vert(VectorLike &res, const Func<Real(const Point<3> &)> &u, const Int &vid) const;
        template <typename VectorLike>
        void interp_edge(VectorLike &res, const Func<Real(const Point<3> &)> &u, const Int &eid) const;
        template <typename VectorLike>
        void interp_face(VectorLike &res, const Func<Real(const Point<3> &)> &u, const Int &fid) const;
        template <typename VectorLike>
        void interp_poly(VectorLike &res, const Func<Real(const Point<3> &)> &u, const Int &pid) const;

        Int get_num_dof_poly(const Int &pid, const Index<3> &i) const;

        Tensor<Real, 2> get_value(const Func<Real(const Point<3> &p)> &u) const;
        Tensor<Real, 2> get_value(const Vector<Real> &u) const;

        Real inner_product(const Vector<Real> &x, const Vector<Real> &y) const;
        Real inner_product(const Tensor<Real, 2> &x, const Tensor<Real, 2> &y) const;
        Real l2_norm(const Vector<Real> &x) const;
        Real l2_norm(const Tensor<Real, 2> &x) const;
        Vector<Real> l2_normalize(const Vector<Real> &x) const;
        Tensor<Real, 2> l2_normalize(const Tensor<Real, 2> &x) const;

        Space &operator=(const Space &space) = delete;
    };
};

#endif
