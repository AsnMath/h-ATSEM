#ifndef ASN_MATH_SPACE_ELLIPTIC_ADAPTER_HPP
#define ASN_MATH_SPACE_ELLIPTIC_ADAPTER_HPP

#include "EllipticAdapter.h"

#include "../Adapter.hpp"

namespace Asn::Math
{
    template <typename Space>
    inline EllipticAdapter<Space>::EllipticAdapter(Space &space) : Adapter<Space>(space) {}

    template <typename Space>
    inline void EllipticAdapter<Space>::set_solution(const Matrix<Real> &solution)
    {
        this->solution = solution;
        this->jump_term = this->build_jump_term(solution);
        return;
    }

    template <typename Space>
    inline void EllipticAdapter<Space>::set_rhs(const Tensor<Real, 2> &rhs)
    {
        this->rhs = rhs;
        return;
    }

    template <typename Space>
    inline Real EllipticAdapter<Space>::h_indicator(const Int &pid) const
    {
        Real interior_residual = 0.0;
        Real boundary_residual = 0.0;

        const Matrix<Real, 3, 3> &inv_jacobi_matrix = this->space.get_poly(pid).inv_jacobi;
        const Matrix<Real, 3, 3> inv_jacobi_matrix_transpose = inv_jacobi_matrix.transpose();

        for (Int l = 0; l < static_cast<Int>(this->space.QUAD.tet.size()); l++)
        {
            Real tmp = 0.0;
            tmp += rhs(pid, l);
            for (const Index<3> &index : this->space.INDEX_POLY_ALL)
            {
                tmp += solution(this->space.get_num_dof_poly(pid, index)) * (inv_jacobi_matrix_transpose * this->space.hess_psi_value(index)[l] * inv_jacobi_matrix).trace();
            }
            interior_residual += (this->space.get_poly(pid).volume * tmp * tmp);
        }

        const Real h_p = this->space.get_poly_diameter(pid) / static_cast<Real>(this->space.MAX_ORDER);
        interior_residual *= (h_p * h_p);

        for (const Int fid : this->space.get_poly(pid).face)
        {
            boundary_residual += jump_term(0, fid);
        }

        return std::sqrt(interior_residual + boundary_residual / 2.0);
    }
};

#endif
