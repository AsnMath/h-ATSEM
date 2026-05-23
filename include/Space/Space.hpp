#ifndef ASN_MATH_SPACE_SPACE_HPP
#define ASN_MATH_SPACE_SPACE_HPP

#include "Space.h"

#include "Index.hpp"

#include "../Quadrature/Quadrature.hpp"
#include "../Mesh/TetMesh.hpp"

namespace Asn::Math
{
    template <Real TOL>
    template <typename Type>
    inline Type Space<TOL>::jacobi(const Int &a, const Int &b, const Int &k, const Type &tu, const Type &td)
    {
        if (k == 0)
        {
            return 1.0;
        }
        if (a == -1 && b == -1)
        {
            if (k == 1)
            {
                return tu;
            }
            return (tu + td) * (tu - td) * jacobi(1, 1, k - 2, tu, td) / 4.0;
        }
        if (a == -1)
        {
            return (k + b) * (tu - td) * jacobi(1, b, k - 1, tu, td) / (2.0 * static_cast<Real>(k));
        }
        if (b == -1)
        {
            return (k + a) * (tu + td) * jacobi(a, 1, k - 1, tu, td) / (2.0 * static_cast<Real>(k));
        }
        Type res = 0.0;
        List<Type> x(k + 1);
        List<Type> y(k + 1);
        x[0] = 1.0;
        y[0] = 1.0;
        for (Int j = 1; j <= k; j++)
        {
            x[j] = x[j - 1] * (tu - td) / 2.0;
            y[j] = y[j - 1] * td;
        }
        for (Int j = 0; j <= k; j++)
        {
            Type tmp = 1.0;
            for (Int i = 1; i <= j; i++)
            {
                tmp = tmp * (k + a + b + i);
                tmp = tmp / i;
            }
            for (Int i = 1; i <= (k - j); i++)
            {
                tmp = tmp * (a + j + i);
                tmp = tmp / i;
            }
            tmp = tmp * x[j] * y[k - j];
            res = res + tmp;
        }
        return res;
    }

    template <Real TOL>
    template <typename Type>
    inline Type Space<TOL>::koornwinder(const Index<3> &i, const Point<3, Type> &x)
    {
        Type res = jacobi<Type>(-1, -1, i(0), 2.0 * x(0) + x(1) + x(2) - 1.0, 1.0 - x(1) - x(2));
        res = res * jacobi<Type>(2 * i(0) - 1, -1, i(1), 2.0 * x(1) + x(2) - 1.0, 1.0 - x(2));
        res = res * jacobi<Type>(2 * i(0) + 2 * i(1) - 1, -1, i(2), 2.0 * x(2) - 1.0, 1.0);
        return res;
    }

    template <Real TOL>
    template <typename Type>
    inline Type Space<TOL>::psi(const Index<3> &i, const Point<3, Type> &x)
    {
        if (i(0) == 0 && i(1) == 0 && i(2) == 0)
        {
            Type res = koornwinder({0, 0, 0}, x) / 8.0;
            res = res - koornwinder({1, 0, 0}, x) / 2.0;
            res = res - koornwinder({0, 1, 0}, x) / 4.0;
            res = res - koornwinder({0, 0, 1}, x) / 8.0;
            return res;
        }
        if (i(0) == 1 && i(1) == 0 && i(2) == 0)
        {
            Type res = koornwinder({0, 0, 0}, x) / 8.0;
            res = res + koornwinder({1, 0, 0}, x) / 2.0;
            res = res - koornwinder({0, 1, 0}, x) / 4.0;
            res = res - koornwinder({0, 0, 1}, x) / 8.0;
            return res;
        }
        if (i(0) == 0 && i(1) == 1 && i(2) == 0)
        {
            Type res = koornwinder({0, 0, 0}, x) / 4.0;
            res = res + koornwinder({0, 1, 0}, x) / 2.0;
            res = res - koornwinder({0, 0, 1}, x) / 4.0;
            return res;
        }
        if (i(0) == 0 && i(1) == 0 && i(2) == 1)
        {
            Type res = koornwinder({0, 0, 0}, x) / 2.0;
            res = res + koornwinder({0, 0, 1}, x) / 2.0;
            return res;
        }

        if (i(0) >= 2 && i(1) == 0 && i(2) == 0)
        {
            Type res = 2.0 * koornwinder({i(0), 0, 0}, x);
            return res;
        }
        if (i(0) == 0 && i(1) >= 2 && i(2) == 0)
        {
            Type res = koornwinder({0, i(1), 0}, x);
            res = res + (i(1) - 1.0) * koornwinder({1, i(1) - 1, 0}, x) / i(1);
            return res;
        }
        if (i(0) == 0 && i(1) == 0 && i(2) >= 2)
        {
            Type res = koornwinder({0, 0, i(2)}, x) / 2.0;
            res = res + (i(2) - 1.0) * koornwinder({0, 1, i(2) - 1}, x) / (2.0 * i(2));
            res = res + (i(2) - 1.0) * koornwinder({1, 0, i(2) - 1}, x) / i(2);
            return res;
        }
        if (i(0) == 1 && i(1) >= 1 && i(2) == 0)
        {
            Type res = koornwinder({0, i(1) + 1, 0}, x);
            res = res - i(1) * koornwinder({1, i(1), 0}, x) / (i(1) + 1);
            return res;
        }
        if (i(0) == 1 && i(1) == 0 && i(2) >= 1)
        {
            Type res = koornwinder({0, 0, i(2) + 1}, x) / 2.0;
            res = res + i(2) * koornwinder({0, 1, i(2)}, x) / (2.0 * (i(2) + 1));
            res = res - i(2) * koornwinder({1, 0, i(2)}, x) / (i(2) + 1);
            return res;
        }
        if (i(0) == 0 && i(1) == 1 && i(2) >= 1)
        {
            Type res = koornwinder({0, 0, i(2) + 1}, x);
            res = res - i(2) * koornwinder({0, 1, i(2)}, x) / (i(2) + 1);
            return res;
        }
        
        if (i(0) == 1 && i(1) >= 1 && i(2) >= 1)
        {
            Type res = koornwinder({0, i(1) + 1, i(2)}, x);
            res = res - i(1) * koornwinder({1, i(1), i(2)}, x) / (i(1) + 1);
            return res;
        }
        if (i(0) == 0 && i(1) >= 2 && i(2) >= 1)
        {
            Type res = koornwinder({0, i(1), i(2)}, x);
            res = res + (i(1) - 1) * koornwinder({1, i(1) - 1, i(2)}, x) / i(1);
            return res;
        }
        if (i(0) >= 2 && i(1) == 0 && i(2) >= 1)
        {
            Type res = 2.0 * koornwinder({i(0), 0, i(2)}, x);
            return res;
        }
        if (i(0) >= 2 && i(1) >= 1 && i(2) == 0)
        {
            Type res = 2.0 * koornwinder({i(0), i(1), 0}, x);
            return res;
        }

        return koornwinder(i, x);
    }

    template <Real TOL>
    inline Point<3> Space<TOL>::grad_psi(const Index<3> &i, const Point<3> &x)
    {
        Point<3, Var> t = x;
        const Var f = psi(i, t);
        Point<3> g = gradient(f, t);
        return g;
    }

    template <Real TOL>
    inline Pair<Point<3>, Matrix<Real, 3, 3>> Space<TOL>::grad_hess_psi(const Index<3> &i, const Point<3> &x)
    {
        Point<3, Var> t = x;
        const Var f = psi<Var>(i, t);
        Point<3> grad;
        Matrix<Real, 3, 3> hess = hessian(f, t, grad);
        return std::make_pair(grad, hess);
    }

    template <Real TOL>
    inline List<Index<2>> Space<TOL>::get_index_list_face(const Int &order)
    {
        List<Index<2>> res;
        for (Int i = 2; i <= order; i++)
        {
            for (Int j = 1; j <= order; j++)
            {
                if ((i + j) > order)
                {
                    break;
                }
                res.emplace_back(i, j);
            }
        }
        return res;
    }

    template <Real TOL>
    inline List<Index<3>> Space<TOL>::get_index_list_poly(const Int &order)
    {
        List<Index<3>> res;
        for (Int i = 2; i <= order; i++)
        {
            for (Int j = 1; j <= order; j++)
            {
                if ((i + j) > order)
                {
                    break;
                }
                for (Int k = 1; k <= order; k++)
                {
                    if ((i + j + k) > order)
                    {
                        break;
                    }
                    res.emplace_back(i, j, k);
                }
            }
        }
        return res;
    }

    template <Real TOL>
    inline List<Index<3>> Space<TOL>::get_index_list_poly_all(const Int &order)
    {
        List<Index<3>> res;
        for (Int i = 0; i <= order; i++)
        {
            for (Int j = 0; j <= order; j++)
            {
                if ((i + j) > order)
                {
                    break;
                }
                for (Int k = 0; k <= order; k++)
                {
                    if ((i + j + k) > order)
                    {
                        break;
                    }
                    res.emplace_back(i, j, k);
                }
            }
        }
        return res;
    }

    template <Real TOL>
    inline SparseMatrix<Real> Space<TOL>::map_matrix(const SparseMatrix<Real> &A, const HashMap<Int, Int> &index, const Int &size)
    {
        const Int N = std::max<Int>(index.size(), size);
        SparseMatrix<Real> res(N, N);
        List<List<Eigen::Triplet<Real>>> triplets(A.outerSize());
        #pragma omp parallel for default(none) shared(triplets, A, index)
        for (Int k = 0; k < A.outerSize(); ++k)
        {
            for (SparseMatrix<Real>::InnerIterator it(A, k); it; ++it)
            {
                if (const Int row = it.row(), col = it.col(); index.contains(row) && index.contains(col))
                {
                    triplets[k].emplace_back(index.at(row), index.at(col), it.value());
                }
            }
        }
        List<Eigen::Triplet<Real>> list;
        Int total_size = 0;
        for (Int i = 0; i < A.outerSize(); i++)
        {
            total_size += triplets[i].size();
        }
        list.reserve(total_size);
        for (Int i = 0; i < A.outerSize(); i++)
        {
            list.insert(list.end(), triplets[i].begin(), triplets[i].end());
        }
        res.setFromTriplets(list.begin(), list.end());
        res.makeCompressed();
        return res;
    }

    template <Real TOL>
    inline Matrix<Real> Space<TOL>::map_matrix_row(const Matrix<Real> &A, const HashMap<Int, Int> &index, const Int &row)
    {
        const Int N = std::max<Int>(index.size(), row);
        Matrix<Real> res = Matrix<Real>::Zero(N, A.cols());
        #pragma omp parallel for default(none) shared(res, A, index)
        for (Int i = 0; i < A.rows(); i++)
        {
            if (index.contains(i))
            {
                res.row(index.at(i)) = A.row(i);
            }
        }
        return res;
    }

    template <Real TOL>
    inline void Space<TOL>::build_dof()
    {
        _num_dof_vert.clear();
        _num_dof_vert.resize(this->num_vert());
        _num_dof_edge.clear();
        _num_dof_edge.resize(this->num_edge());
        _num_dof_face.clear();
        _num_dof_face.resize(this->num_face());
        _num_dof_poly.clear();
        _num_dof_poly.resize(this->num_poly());
        _num_dof_poly_all.clear();
        _num_dof_poly_all.resize(this->num_poly());
        #pragma omp parallel for default(none)
        for (Int fid = 0; fid < this->num_face(); fid++)
        {
            _num_dof_face[fid].init(MAX_ORDER, -1);
        }
        #pragma omp parallel for default(none)
        for (Int pid = 0; pid < this->num_poly(); pid++)
        {
            _num_dof_poly[pid].init(MAX_ORDER, -1);
        }

        _num_dof = 0;

        for (Int vid = 0; vid < this->num_vert(); vid++)
        {
            _num_dof_vert[vid] = _num_dof;
            _num_dof++;
        }

        for (Int eid = 0; eid < this->num_edge(); eid++)
        {
            _num_dof_edge[eid].resize(MAX_ORDER + 1);
            for (Int order = 2; order <= MAX_ORDER; order++)
            {
                _num_dof_edge[eid][order] = _num_dof;
                _num_dof++;
            }
        }
        for (Int fid = 0; fid < this->num_face(); fid++)
        {
            for (const Index<2> &index : INDEX_FACE)
            {
                _num_dof_face[fid](index) = _num_dof;
                _num_dof++;
            }
        }
        for (Int pid = 0; pid < this->num_poly(); pid++)
        {
            for (const Index<3> &index : INDEX_POLY)
            {
                _num_dof_poly[pid](index) = _num_dof;
                _num_dof++;
            }
        }

        #pragma omp parallel for default(none)
        for (Int pid = 0; pid < this->num_poly(); pid++)
        {
            _num_dof_poly_all[pid].resize(INDEX_POLY_ALL.size());
            for (Int it = 0; it < static_cast<Int>(INDEX_POLY_ALL.size()); it++)
            {
                _num_dof_poly_all[pid][it] = get_num_dof_poly(pid, INDEX_POLY_ALL[it]);
            }
        }
        return;
    }

    template <Real TOL>
    inline void Space<TOL>::build_precomputation()
    {
        const Array<Point<3>, 4> p = {
                    Point<3>({0.0, 0.0, 0.0}),
                    Point<3>({1.0, 0.0, 0.0}),
                    Point<3>({0.0, 1.0, 0.0}),
                    Point<3>({0.0, 0.0, 1.0})
                };
        const Int N_POINT_1 = static_cast<Int>(QUAD.line.size());
        const Int N_POINT_2 = static_cast<Int>(QUAD.tri.size());
        const Int N_POINT_3 = static_cast<Int>(QUAD.tet.size());

        _psi_value.init(MAX_ORDER, List<Real>());
        _grad_psi_value.init(MAX_ORDER, List<Point<3>>());
        _grad_psi_value_face.init(MAX_ORDER, Array<List<Point<3>>, 4>());
        _hess_psi_value.init(MAX_ORDER, List<Matrix<Real, 3, 3>>());

        #pragma omp parallel for default(none) shared(p, N_POINT_2, N_POINT_3)
        for (Int i = 0; i < static_cast<Int>(INDEX_POLY_ALL.size()); i++)
        {
            const Index<3> &index = INDEX_POLY_ALL[i];
            _psi_value(index).resize(N_POINT_3);
            _grad_psi_value(index).resize(N_POINT_3);
            _hess_psi_value(index).resize(N_POINT_3);
            for (Int l = 0; l < N_POINT_3; l++)
            {
                const Point<3> &x = QUAD.tet[l].first;
                _psi_value(index)[l] = psi<Real>(index, x);
                std::tie(_grad_psi_value(index)[l], _hess_psi_value(index)[l]) = grad_hess_psi(index, x);
            }
            _grad_psi_value_face(index)[0].resize(N_POINT_2);
            _grad_psi_value_face(index)[1].resize(N_POINT_2);
            _grad_psi_value_face(index)[2].resize(N_POINT_2);
            _grad_psi_value_face(index)[3].resize(N_POINT_2);
            for (Int l = 0; l < N_POINT_2; l++)
            {
                const Point<2> &x = QUAD.tri[l].first;
                _grad_psi_value_face(index)[0][l] = grad_psi(index, p[1] + x(0) * (p[2] - p[1]) + x(1) * (p[3] - p[1]));
                _grad_psi_value_face(index)[1][l] = grad_psi(index, p[0] + x(0) * (p[2] - p[0]) + x(1) * (p[3] - p[0]));
                _grad_psi_value_face(index)[2][l] = grad_psi(index, p[0] + x(0) * (p[1] - p[0]) + x(1) * (p[3] - p[0]));
                _grad_psi_value_face(index)[3][l] = grad_psi(index, p[0] + x(0) * (p[1] - p[0]) + x(1) * (p[2] - p[0]));
            }
        }
        _inner_product_value.resize(INDEX_POLY_ALL.size());
        #pragma omp parallel for default(none)
        for (Int i = 0; i < static_cast<Int>(INDEX_POLY_ALL.size()); i++)
        {
            _inner_product_value[i].resize(INDEX_POLY_ALL.size(), 0.0);
            const Index<3> &index_1 = INDEX_POLY_ALL[i];
            for (Int j = i; j < static_cast<Int>(INDEX_POLY_ALL.size()); j++)
            {
                const Index<3> &index_2 = INDEX_POLY_ALL[j];
                for (Int l = 0; l < static_cast<Int>(QUAD.tet.size()); l++)
                {
                    const Real &w = QUAD.tet[l].second;
                    _inner_product_value[i][j] += (w * psi_value(index_1)[l] * psi_value(index_2)[l]);
                }
            }
        }
        #pragma omp parallel for default(none)
        for (Int i = 0; i < static_cast<Int>(INDEX_POLY_ALL.size()); i++)
        {
            for (Int j = 0; j < i; j++)
            {
                _inner_product_value[i][j] = _inner_product_value[j][i];
            }
        }

        jacobi_edge.resize(MAX_ORDER + 1);
        jacobi_face_1.resize(MAX_ORDER + 1);
        jacobi_face_2.resize(MAX_ORDER + 1);
        #pragma omp parallel for default(none) shared(N_POINT_1, N_POINT_2)
        for (Int k = 0; k <= MAX_ORDER; k++)
        {
            jacobi_edge[k].resize(N_POINT_1);
            jacobi_face_1[k].resize(N_POINT_2);
            jacobi_face_2[k].resize(N_POINT_2);
        }
        #pragma omp parallel for default(none) shared(N_POINT_1)
        for (Int l = 0; l < N_POINT_1; l++)
        {
            const Real &x = QUAD.line[l].first(0);
            for (Int k = 0; k <= MAX_ORDER; k++)
            {
                jacobi_edge[k][l] = jacobi<Real>(1, 1, k, 2.0 * x - 1.0, 1.0);
            }
        }
        #pragma omp parallel for default(none) shared(N_POINT_2)
        for (Int l = 0; l < N_POINT_2; l++)
        {
            const Real &x = QUAD.tri[l].first(0);
            const Real &y = QUAD.tri[l].first(1);
            for (Int k = 0; k <= MAX_ORDER; k++)
            {
                jacobi_face_1[k][l] = jacobi<Real>(1, 1, k, 2.0 * x + y - 1.0, 1.0 - y);
                jacobi_face_2[k][l] = jacobi<Real>(1, 1, k, 2.0 * y - 1.0, 1.0);
            }
        }
        return;
    }

    template <Real TOL>
    inline void Space<TOL>::build_precomputation_element()
    {
        const Int N_POINT_3 = static_cast<Int>(QUAD.tet.size());

        _trans_grad_psi_value.resize(this->num_poly());

        #pragma omp parallel for default(none) shared(N_POINT_3)
        for (Int pid = 0; pid < this->num_poly(); pid++)
        {
            const Matrix<Real, 3, 3> inv_jacobi_matrix_transpose = this->get_poly(pid).inv_jacobi.transpose();
            _trans_grad_psi_value[pid].init(MAX_ORDER, List<Point<3>>(QUAD.tet.size()));
            for (const Index<3> &index : INDEX_POLY_ALL)
            {
                for (Int l = 0; l < N_POINT_3; l++)
                {
                    _trans_grad_psi_value[pid](index)[l] = inv_jacobi_matrix_transpose * _grad_psi_value(index)[l];
                }
            }
        }
        return;
    }

    template <Real TOL>
    inline void Space<TOL>::build_boundary()
    {
        _boundary_vert.clear();
        _boundary_edge.clear();
        _boundary_face.clear();
        _boundary_dof.clear();

        for (Int fid = 0; fid < this->num_face(); fid++)
        {
            if (this->get_face(fid).poly.size() == 1)
            {
                _boundary_face.push_back(fid);
                const List<Int> edge = convert<HashSet<Int>, List<Int>>(this->get_face(fid).edge);
                const Array<Int, BaseMesh::VERTS_PER_FACE> vert = this->get_face(fid).vert;
                for (const Int &item : edge)
                {
                    _boundary_edge.push_back(item);
                }
                for (const Int &item : vert)
                {
                    _boundary_vert.push_back(item);
                }
            }
        }

        std::ranges::sort(_boundary_face);
        _boundary_face.erase(std::ranges::unique(_boundary_face).begin(), _boundary_face.end());
        std::ranges::sort(_boundary_edge);
        _boundary_edge.erase(std::ranges::unique(_boundary_edge).begin(), _boundary_edge.end());
        std::ranges::sort(_boundary_vert);
        _boundary_vert.erase(std::ranges::unique(_boundary_vert).begin(), _boundary_vert.end());

        for (Int i = 0; i < static_cast<Int>(boundary_vert.size()); i++)
        {
            _boundary_dof.insert(num_dof_vert[boundary_vert[i]]);
        }

        for (Int i = 0; i < static_cast<Int>(boundary_edge.size()); i++)
        {
            const Int eid = boundary_edge[i];
            for (Int j = 2; j <= MAX_ORDER; j++)
            {
                _boundary_dof.insert(num_dof_edge[eid][j]);
            }
        }

        for (Int i = 0; i < static_cast<Int>(boundary_face.size()); i++)
        {
            const Int fid = boundary_face[i];
            for (const Index<2> &index : INDEX_FACE)
            {
                _boundary_dof.insert(num_dof_face[fid](index));
            }
        }

        return;
    }

    template <Real TOL>
    inline Space<TOL>::Space(const String &mesh_file, const Int &order, const Int &quad_order)
        : BaseMesh(),
          MAX_ORDER(order), QUAD_ORDER(quad_order), QUAD(QUAD_ORDER),
          INDEX_FACE(get_index_list_face(order)), INDEX_POLY(get_index_list_poly(order)),
          INDEX_POLY_ALL(get_index_list_poly_all(order)),
          num_dof(_num_dof), num_dof_vert(_num_dof_vert), num_dof_edge(_num_dof_edge),
          num_dof_face(_num_dof_face), num_dof_poly(_num_dof_poly),
          num_dof_poly_all(_num_dof_poly_all),
          inner_product_value(_inner_product_value), psi_value(_psi_value),
          grad_psi_value(_grad_psi_value), grad_psi_value_face(_grad_psi_value_face),
          hess_psi_value(_hess_psi_value), trans_grad_psi_value(_trans_grad_psi_value),
          boundary_vert(_boundary_vert), boundary_edge(_boundary_edge),
          boundary_face(_boundary_face), boundary_dof(_boundary_dof),
          stiff_matrix(_stiff_matrix), mass_matrix(_mass_matrix)
    {
        assert(quad_order >= std::max<Int>(2 * (order - 1), 1));

        this->load(mesh_file);
        this->init_h_adaptive();
        build_precomputation();

        init();
        return;
    }

    template <Real TOL>
    inline Space<TOL>::Space(const Space &space)
        : BaseMesh(space),
          MAX_ORDER(space.MAX_ORDER), QUAD_ORDER(space.QUAD_ORDER), QUAD(space.QUAD),
          INDEX_FACE(space.INDEX_FACE), INDEX_POLY(space.INDEX_POLY),
          INDEX_POLY_ALL(space.INDEX_POLY_ALL),
          num_dof(_num_dof), num_dof_vert(_num_dof_vert), num_dof_edge(_num_dof_edge),
          num_dof_face(_num_dof_face), num_dof_poly(_num_dof_poly),
          num_dof_poly_all(_num_dof_poly_all),
          inner_product_value(_inner_product_value), psi_value(_psi_value),
          grad_psi_value(_grad_psi_value), grad_psi_value_face(_grad_psi_value_face),
          hess_psi_value(_hess_psi_value), trans_grad_psi_value(_trans_grad_psi_value),
          boundary_vert(_boundary_vert), boundary_edge(_boundary_edge),
          boundary_face(_boundary_face), boundary_dof(_boundary_dof),
          stiff_matrix(_stiff_matrix), mass_matrix(_mass_matrix)
    {
        _num_dof = space.num_dof;
        _num_dof_vert = space.num_dof_vert;
        _num_dof_edge = space.num_dof_edge;
        _num_dof_face = space.num_dof_face;
        _num_dof_poly = space.num_dof_poly;
        _num_dof_poly_all = space.num_dof_poly_all;
        _inner_product_value = space.inner_product_value;
        _psi_value = space.psi_value;
        _grad_psi_value = space.grad_psi_value;
        _grad_psi_value_face = space.grad_psi_value_face;
        _hess_psi_value = space.hess_psi_value;
        _trans_grad_psi_value = space.trans_grad_psi_value;
        _boundary_vert = space.boundary_vert;
        _boundary_edge = space.boundary_edge;
        _boundary_face = space.boundary_face;
        _boundary_dof = space.boundary_dof;
    }

    template <Real TOL>
    inline void Space<TOL>::init()
    {
        this->update_attributes();
        build_dof();
        build_boundary();
        build_precomputation_element();
        return;
    }

    template <Real TOL>
    inline Real Space<TOL>::build_element_stiff_matrix(const Int &pid, const Index<3> &i, const Index<3> &j) const
    {
        Real res = 0.0;
        for (Int l = 0; l < static_cast<Int>(QUAD.tet.size()); l++)
        {
            const Point<3> &g1 = trans_grad_psi_value[pid](i)[l];
            const Point<3> &g2 = trans_grad_psi_value[pid](j)[l];
            res += QUAD.tet[l].second * g1.dot(g2);
        }
        return this->poly[pid].volume * res;
    }

    template <Real TOL>
    inline void Space<TOL>::build_stiff_matrix()
    {
        List<List<Eigen::Triplet<Real>>> list(this->poly.size());
        #pragma omp parallel for default(none) shared(list)
        for (Int pid = 0; pid < static_cast<Int>(this->poly.size()); pid++)
        {
            for (Int it = 0; it < static_cast<Int>(INDEX_POLY_ALL.size()); it++)
            {
                const Int dof_1 = num_dof_poly_all[pid][it];
                list[pid].emplace_back(dof_1, dof_1, build_element_stiff_matrix(pid, INDEX_POLY_ALL[it], INDEX_POLY_ALL[it]));
                for (Int jt = it + 1; jt < static_cast<Int>(INDEX_POLY_ALL.size()); jt++)
                {
                    if (const Int dof_2 = num_dof_poly_all[pid][jt]; dof_1 > dof_2)
                    {
                        list[pid].emplace_back(dof_1, dof_2, build_element_stiff_matrix(pid, INDEX_POLY_ALL[it], INDEX_POLY_ALL[jt]));
                    }
                    else
                    {
                        list[pid].emplace_back(dof_2, dof_1, build_element_stiff_matrix(pid, INDEX_POLY_ALL[it], INDEX_POLY_ALL[jt]));
                    }
                }
            }
        }
        List<Eigen::Triplet<Real>> list_all;
        Int size = 0;
        for (Int pid = 0; pid < static_cast<Int>(this->poly.size()); pid++)
        {
            size += list[pid].size();
        }
        list_all.reserve(size);
        for (Int pid = 0; pid < static_cast<Int>(this->poly.size()); pid++)
        {
            list_all.insert(list_all.end(), list[pid].begin(), list[pid].end());
        }
        SparseMatrix<Real> A(num_dof, num_dof);
        A.setFromTriplets(list_all.begin(), list_all.end());
        A.prune(static_cast<Real>(0.0));
        _stiff_matrix = A.selfadjointView<Eigen::Lower>();
        _stiff_matrix.makeCompressed();
        return;
    }

    template <Real TOL>
    inline void Space<TOL>::build_mass_matrix()
    {
        List<List<Eigen::Triplet<Real>>> list(this->poly.size());
        #pragma omp parallel for default(none) shared(list)
        for (Int pid = 0; pid < static_cast<Int>(this->poly.size()); pid++)
        {
            for (Int it = 0; it < static_cast<Int>(INDEX_POLY_ALL.size()); it++)
            {
                const Real volume = this->poly[pid].volume;
                const Int dof_1 = num_dof_poly_all[pid][it];
                list[pid].emplace_back(dof_1, dof_1, volume * inner_product_value[it][it]);
                for (Int jt = it + 1; jt < static_cast<Int>(INDEX_POLY_ALL.size()); jt++)
                {
                    if (const Int dof_2 = num_dof_poly_all[pid][jt]; dof_1 > dof_2)
                    {
                        list[pid].emplace_back(dof_1, dof_2, volume * inner_product_value[it][jt]);
                    }
                    else
                    {
                        list[pid].emplace_back(dof_2, dof_1, volume * inner_product_value[it][jt]);
                    }
                }
            }
        }
        List<Eigen::Triplet<Real>> list_all;
        Int size = 0;
        for (Int pid = 0; pid < static_cast<Int>(this->poly.size()); pid++)
        {
            size += list[pid].size();
        }
        list_all.reserve(size);
        for (Int pid = 0; pid < static_cast<Int>(this->poly.size()); pid++)
        {
            list_all.insert(list_all.end(), list[pid].begin(), list[pid].end());
        }
        SparseMatrix<Real> A(num_dof, num_dof);
        A.setFromTriplets(list_all.begin(), list_all.end());
        A.prune(static_cast<Real>(0.0));
        _mass_matrix = A.selfadjointView<Eigen::Lower>();
        _mass_matrix.makeCompressed();
        return;
    }

    template <Real TOL>
    template <typename VectorLike>
    inline void Space<TOL>::interp_vert(VectorLike &res, const Func<Real(const Point<3> &)> &u, const Int &vid) const
    {
        const Point<3> &p = this->vert[vid].coord;
        res(num_dof_vert[vid]) = u(p);
        return;
    }

    template <Real TOL>
    template <typename VectorLike>
    inline void Space<TOL>::interp_edge(VectorLike &res, const Func<Real(const Point<3> &)> &u, const Int &eid) const
    {
        const Array<Int, BaseMesh::VERTS_PER_EDGE> v = this->get_edge(eid).vert;
        const Point<3> &ps = this->vert[v[0]].coord;
        const Point<3> &pe = this->vert[v[1]].coord;

        List<Real> coef(MAX_ORDER + 1, 0.0);
        coef[0] = res(num_dof_vert[v[0]]);
        coef[1] = res(num_dof_vert[v[1]]);

        List<Real> func_value(QUAD.line.size());
        for (Int l = 0; l < static_cast<Int>(QUAD.line.size()); l++)
        {
            const Real &x = QUAD.line[l].first(0);
            func_value[l] = u((1.0 - x) * ps + x * pe) - coef[0] * (1 - x) - coef[1] * x;
        }

        for (Int k = 2; k <= MAX_ORDER; k++)
        {
            for (Int l = 0; l < static_cast<Int>(QUAD.line.size()); l++)
            {
                const Real &w = QUAD.line[l].second;
                coef[k] += w * func_value[l] * jacobi_edge[k - 2][l];
            }
            coef[k] *= -(k * (2.0 * k - 1.0)) / (2.0 * (k - 1.0));
            res[num_dof_edge[eid][k]] = coef[k];
        }
        return;
    }

    template <Real TOL>
    template <typename VectorLike>
    inline void Space<TOL>::interp_face(VectorLike &res, const Func<Real(const Point<3> &)> &u, const Int &fid) const
    {
        const Array<Int, BaseMesh::VERTS_PER_FACE> vid = this->get_face(fid).vert;
        const Real c_0 = res(vid[0]);
        const Real c_1 = res(vid[1]);
        const Real c_2 = res(vid[2]);
        const Int eid_01 = this->find_edge({vid[0], vid[1]});
        const Int eid_02 = this->find_edge({vid[0], vid[2]});
        const Int eid_12 = this->find_edge({vid[1], vid[2]});
        List<Real> c_01(MAX_ORDER + 1, 0.0);
        List<Real> c_02(MAX_ORDER + 1, 0.0);
        List<Real> c_12(MAX_ORDER + 1, 0.0);
        for (Int i = 2; i <= MAX_ORDER; i++)
        {
            c_01[i] = res(num_dof_edge[eid_01][i]);
            c_02[i] = res(num_dof_edge[eid_02][i]);
            c_12[i] = res(num_dof_edge[eid_12][i]);
        }
        const Point<3> &p_0 = this->vert[vid[0]].coord;
        const Point<3> &p_1 = this->vert[vid[1]].coord;
        const Point<3> &p_2 = this->vert[vid[2]].coord;

        List<Real> func_value(QUAD.tri.size());
        for (Int l = 0; l < static_cast<Int>(QUAD.tri.size()); l++)
        {
            const Point<2> &x = QUAD.tri[l].first;
            func_value[l] = u(p_0 + (p_1 - p_0) * x(0) + (p_2 - p_0) * x(1)) - (c_0 * (1.0 - x(0) - x(1)) + c_1 * x(0) + c_2 * x(1));

            Real tmp_01 = 0.0;
            Real tmp_02 = 0.0;
            Real tmp_12 = 0.0;

            for (Int m = 2; m <= MAX_ORDER; m++)
            {
                tmp_01 += c_01[m] * jacobi_face_1[m - 2][l];
                tmp_02 += c_02[m] * jacobi_face_2[m - 2][l];
                tmp_12 += c_12[m] * jacobi_face_2[m - 2][l];
            }

            const Real r = 1.0 - x(0) - x(1);
            tmp_01 *= (2.0 * x(0) * r);
            tmp_02 *= (2.0 * x(1) * r);
            tmp_12 *= (2.0 * x(0) * x(1));

            func_value[l] += tmp_01 + tmp_02 + tmp_12;
        }

        for (const Index<2> &index : INDEX_FACE)
        {
            const Int &i = index(0);
            const Int &j = index(1);
            Real &c = res(num_dof_face[fid](index));
            c = 0.0;
            for (Int l = 0; l < static_cast<Int>(QUAD.tri.size()); l++)
            {
                const Point<2> &p = QUAD.tri[l].first;
                const Real &w = QUAD.tri[l].second;
                const Real x = p(0);
                const Real y = p(1);
                const Real r = 1.0 - x - y;

                Real tmp = func_value[l];
                tmp *= jacobi_face_1[i - 2][l] * jacobi<Real>(2 * i - 1, 1, j - 1, 2.0 * y - 1.0, 1.0);
                c += w * tmp;
            }
            c *= (-i * (2.0 * i - 1.0) * (2.0 * i + 2.0 * j - 1.0)) / (4.0 * (i - 1.0));
        }
        return;
    }

    template <Real TOL>
    template <typename VectorLike>
    inline void Space<TOL>::interp_poly(VectorLike &res, const Func<Real(const Point<3> &)> &u, const Int &pid) const
    {
        const Point<3> &p_0 = this->vert[this->get_poly(pid).vert[0]].coord;
        const Matrix<Real, 3, 3> &jacobi_matrix = this->poly[pid].jacobi;
        List<Real> func_value(QUAD.tet.size());
        for (Int l = 0; l < static_cast<Int>(QUAD.tet.size()); l++)
        {
            const Point<3> &x = QUAD.tet[l].first;
            func_value[l] = u(p_0 + jacobi_matrix * x);
            for (const Index<3> &index_2 : INDEX_POLY_ALL)
            {
                if (index_2(0) < 2 || index_2(1) < 1 || index_2(2) < 1)
                {
                    func_value[l] -= (res(get_num_dof_poly(pid, index_2)) * psi_value(index_2)[l]);
                }
            }
        }

        for (const Index<3> &index : INDEX_POLY)
        {
            Real &c = res(num_dof_poly[pid](index));
            const Int &i = index(0);
            const Int &j = index(1);
            const Int &k = index(2);
            c = 0.0;
            for (Int l = 0; l < static_cast<Int>(QUAD.tet.size()); l++)
            {
                const Point<3> &p = QUAD.tet[l].first;
                const Real &w = QUAD.tet[l].second;
                const Real x = p(0);
                const Real y = p(1);
                const Real z = p(2);
                Real tmp = func_value[l];
                tmp *= jacobi<Real>(1, 1, i - 2, 2.0 * x + y + z - 1.0, 1.0 - y - z);
                tmp *= jacobi<Real>(2 * i - 1, 1, j - 1, 2.0 * y + z - 1.0, 1.0 - z);
                tmp *= jacobi<Real>(2 * i + 2 * j - 1, 1, k - 1, 2.0 * z - 1.0, 1.0);
                c += w * tmp;
            }
            c *= -(i * (2.0 * i - 1.0) * (2.0 * i + 2.0 * j - 1.0) * (2.0 * i + 2.0 * j + 2.0 * k - 1.0)) / (6.0 * (i - 1.0));
        }
        return;
    }

    template <Real TOL>
    inline Int Space<TOL>::get_num_dof_poly(const Int &pid, const Index<3> &i) const
    {
        const Array<Int, BaseMesh::VERTS_PER_POLY> vid = this->get_poly(pid).vert;

        if (i(0) == 0 && i(1) == 0 && i(2) == 0)
        {
            return _num_dof_vert[vid[0]];
        }
        if (i(0) == 1 && i(1) == 0 && i(2) == 0)
        {
            return _num_dof_vert[vid[1]];
        }
        if (i(0) == 0 && i(1) == 1 && i(2) == 0)
        {
            return _num_dof_vert[vid[2]];
        }
        if (i(0) == 0 && i(1) == 0 && i(2) == 1)
        {
            return _num_dof_vert[vid[3]];
        }

        if (i(0) >= 2 && i(1) == 0 && i(2) == 0)
        {
            const Int eid = this->find_edge({vid[0], vid[1]});
            assert(eid >= 0);
            return _num_dof_edge[eid][i(0) + i(1) + i(2)];
        }
        if (i(0) == 0 && i(1) >= 2 && i(2) == 0)
        {
            const Int eid = this->find_edge({vid[0], vid[2]});
            assert(eid >= 0);
            return _num_dof_edge[eid][i(0) + i(1) + i(2)];
        }
        if (i(0) == 0 && i(1) == 0 && i(2) >= 2)
        {
            const Int eid = this->find_edge({vid[0], vid[3]});
            assert(eid >= 0);
            return _num_dof_edge[eid][i(0) + i(1) + i(2)];
        }
        if (i(0) == 1 && i(1) >= 1 && i(2) == 0)
        {
            const Int eid = this->find_edge({vid[1], vid[2]});
            assert(eid >= 0);
            return _num_dof_edge[eid][i(0) + i(1) + i(2)];
        }
        if (i(0) == 1 && i(1) == 0 && i(2) >= 1)
        {
            const Int eid = this->find_edge({vid[1], vid[3]});
            assert(eid >= 0);
            return _num_dof_edge[eid][i(0) + i(1) + i(2)];
        }
        if (i(0) == 0 && i(1) == 1 && i(2) >= 1)
        {
            const Int eid = this->find_edge({vid[2], vid[3]});
            assert(eid >= 0);
            return _num_dof_edge[eid][i(0) + i(1) + i(2)];
        }

        if (i(0) == 1 && i(1) >= 1 && i(2) >= 1)
        {
            const Int fid = this->find_face({vid[1], vid[2], vid[3]});
            assert(fid >= 0);
            return _num_dof_face[fid](Index<2>({i(1) + 1, i(2)}));
        }
        if (i(0) == 0 && i(1) >= 2 && i(2) >= 1)
        {
            const Int fid = this->find_face({vid[0], vid[2], vid[3]});
            assert(fid >= 0);
            return _num_dof_face[fid](Index<2>({i(1), i(2)}));
        }
        if (i(0) >= 2 && i(1) == 0 && i(2) >= 1)
        {
            const Int fid = this->find_face({vid[0], vid[1], vid[3]});
            assert(fid >= 0);
            return _num_dof_face[fid](Index<2>({i(0), i(2)}));
        }
        if (i(0) >= 2 && i(1) >= 1 && i(2) == 0)
        {
            const Int fid = this->find_face({vid[0], vid[1], vid[2]});
            assert(fid >= 0);
            return _num_dof_face[fid](Index<2>({i(0), i(1)}));
        }

        return _num_dof_poly[pid](i);
    }

    template <Real TOL>
    inline Tensor<Real, 2> Space<TOL>::get_value(const Func<Real(const Point<3> &p)> &u) const
    {
        Tensor<Real, 2> res(this->poly.size(), QUAD.tet.size());
        res.setZero();
        #pragma omp parallel for default(none) shared(u, res)
        for (Int pid = 0; pid < static_cast<Int>(this->poly.size()); pid++)
        {
            const Matrix<Real, 3, 3> &jacobi_matrix = this->poly[pid].jacobi;
            const Point<3> &p_0 = this->vert[this->get_poly(pid).vert[0]].coord;
            for (Int l = 0; l < static_cast<Int>(QUAD.tet.size()); l++)
            {
                res(pid, l) = u(p_0 + jacobi_matrix * QUAD.tet[l].first);
            }
        }
        return res;
    }

    template <Real TOL>
    inline Tensor<Real, 2> Space<TOL>::get_value(const Vector<Real> &u) const
    {
        Tensor<Real, 2> res(this->poly.size(), QUAD.tet.size());
        res.setZero();
        #pragma omp parallel for default(none) shared(u, res)
        for (Int pid = 0; pid < static_cast<Int>(this->poly.size()); pid++)
        {
            for (Int it = 0; it < static_cast<Int>(INDEX_POLY_ALL.size()); it++)
            {
                const Index<3> &index = INDEX_POLY_ALL[it];
                const Int dof = num_dof_poly_all[pid][it];
                for (Int l = 0; l < static_cast<Int>(QUAD.tet.size()); l++)
                {
                    res(pid, l) += (u(dof) * psi_value(index)[l]);
                }
            }
        }
        return res;
    }

    template <Real TOL>
    inline Real Space<TOL>::inner_product(const Vector<Real> &x, const Vector<Real> &y) const
    {
        Real res = 0.0;
        List<Real> list_all(this->poly.size(), 0.0);
        #pragma omp parallel for default(none) shared(list_all, x, y)
        for (Int pid = 0; pid < static_cast<Int>(this->poly.size()); pid++)
        {
            for (Int it = 0; it < static_cast<Int>(INDEX_POLY_ALL.size()); it++)
            {
                const Int &dof_1 = num_dof_poly_all[pid][it];
                for (Int jt = 0; jt < static_cast<Int>(INDEX_POLY_ALL.size()); jt++)
                {
                    const Int &dof_2 = num_dof_poly_all[pid][jt];
                    list_all[pid] += x(dof_1) * y(dof_2) * _inner_product_value[it][jt];
                }
            }
            list_all[pid] *= this->poly[pid].volume;
        }
        for (const Real item : list_all)
        {
            res += item;
        }
        return res;
    }

    template <Real TOL>
    inline Real Space<TOL>::inner_product(const Tensor<Real, 2> &x, const Tensor<Real, 2> &y) const
    {
        Real res = 0.0;
        List<Real> list_all(this->poly.size(), 0.0);
        #pragma omp parallel for default(none) shared(list_all, x, y)
        for (Int pid = 0; pid < static_cast<Int>(this->poly.size()); pid++)
        {
            for (Int l = 0; l < static_cast<Int>(QUAD.tet.size()); l++)
            {
                list_all[pid] += QUAD.tet[l].second * x(pid, l) * y(pid, l);
            }
            list_all[pid] *= this->poly[pid].volume;
        }
        for (const Real item : list_all)
        {
            res += item;
        }
        return res;
    }

    template <Real TOL>
    inline Real Space<TOL>::l2_norm(const Vector<Real> &x) const
    {
        const Real norm = std::sqrt(inner_product(x, x));
        return norm;
    }

    template <Real TOL>
    inline Real Space<TOL>::l2_norm(const Tensor<Real, 2> &x) const
    {
        const Real norm = std::sqrt(inner_product(x, x));
        return norm;
    }

    template <Real TOL>
    inline Vector<Real> Space<TOL>::l2_normalize(const Vector<Real> &x) const
    {
        return x / std::sqrt(inner_product(x, x));
    }

    template <Real TOL>
    inline Tensor<Real, 2> Space<TOL>::l2_normalize(const Tensor<Real, 2> &x) const
    {
        return x / std::sqrt(inner_product(x, x));
    }
};

#endif
