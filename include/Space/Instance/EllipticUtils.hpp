#ifndef ASN_MATH_SPACE_ELLIPTIC_UTILS_HPP
#define ASN_MATH_SPACE_ELLIPTIC_UTILS_HPP

#include "EllipticUtils.h"

#include "../Utils.hpp"

namespace Asn::Math
{
    template <typename Space>
    inline EllipticUtils<Space>::EllipticUtils(Space &space) : Utils<Space>(space) {}

    template <typename Space>
    inline void EllipticUtils<Space>::init()
    {
        Utils<Space>::init();
        this->space.build_stiff_matrix();
        return;
    }

    template <typename Space>
    inline Vector<Real> EllipticUtils<Space>::get_boundary(const Func<Real(const Point<Space::MESH_DIM> &)> &u)
    const
    {
        Vector<Real> res = Vector<Real>::Zero(this->space.num_dof);
        #pragma omp parallel for default(none) shared(u, res)
        for (Int i = 0; i < static_cast<Int>(this->space.boundary_vert.size()); i++)
        {
            const Int vid = this->space.boundary_vert[i];
            res(this->space.num_dof_vert[vid]) = u(this->space.get_vert(vid).coord);
        }
        #pragma omp parallel for default(none) shared(u, res)
        for (Int i = 0; i < static_cast<Int>(this->space.boundary_edge.size()); i++)
        {
            this->space.interp_edge(res, u, this->space.boundary_edge[i]);
        }

        #pragma omp parallel for default(none) shared(u, res)
        for (Int i = 0; i < static_cast<Int>(this->space.boundary_face.size()); i++)
        {
            this->space.interp_face(res, u, this->space.boundary_face[i]);
        }
        return res;
    }

    template <typename Space>
    inline Vector<Real> EllipticUtils<Space>::get_rhs(const Tensor<Real, 2> &f) const
    {
        List<List<Pair<Int, Real>>> list(this->space.num_poly());
        #pragma omp parallel for default(none) shared(list, f)
        for (Int pid = 0; pid < this->space.num_poly(); pid++)
        {
            const Real volume = this->space.get_poly(pid).volume;
            for (Int it = 0; it < static_cast<Int>(this->space.INDEX_POLY_ALL.size()); it++)
            {
                const Int dof = this->space.num_dof_poly_all[pid][it];
                if (this->space.boundary_dof.contains(dof))
                {
                    continue;
                }
                Real tmp = 0.0;
                for (Int l = 0; l < static_cast<Int>(this->space.QUAD.tet.size()); l++)
                {
                    tmp += this->space.QUAD.tet[l].second * f(pid, l) * this->space.psi_value(this->space.INDEX_POLY_ALL[it])[l];
                }
                list[pid].emplace_back(dof, volume * tmp);
            }
        }
        Vector<Real> rhs = Vector<Real>::Zero(this->space.num_dof);
        for (Int pid = 0; pid < this->space.num_poly(); pid++)
        {
            for (const auto &[row, data] : list[pid])
            {
                rhs(row) += data;
            }
        }
        return rhs;
    }

    template <typename Space>
    inline void EllipticUtils<Space>::apply_boundary(SparseMatrix<Real> &A, Vector<Real> &u) const
    {
        #pragma omp parallel for default(none) shared(A, u)
        for (Int row = 0; row < this->space.num_dof; row++)
        {
            SparseMatrix<Real>::InnerIterator it(A, row);
            if (this->space.boundary_dof.contains(row))
            {
                for (; it && it.col() < row; ++it)
                {
                    it.valueRef() = 0.0;
                }
                if (it && it.col() == row)
                {
                    it.valueRef() = 1.0;
                    ++it;
                }
                for (; it; ++it)
                {
                    it.valueRef() = 0.0;
                }
            }
            else
            {
                for (; it; ++it)
                {
                    if (const Int col = it.col(); this->space.boundary_dof.contains(it.col()))
                    {
                        u(row) -= (it.value() * u(col));
                        it.valueRef() = 0.0;
                    }
                }
            }
        }
        A.prune(0.0);
        return;
    }
};

#endif
