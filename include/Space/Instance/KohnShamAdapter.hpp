#ifndef ASN_MATH_SPACE_KOHN_SHAM_ADAPTER_HPP
#define ASN_MATH_SPACE_KOHN_SHAM_ADAPTER_HPP

#include <pstl/parallel_backend_utils.h>

#include "KohnShamAdapter.h"

#include "../Adapter.hpp"

namespace Asn::Math
{
    template <typename Space>
    inline KohnShamAdapter<Space>::KohnShamAdapter(Space &space) : Adapter<Space>(space) {}

    template <typename Space>
    inline void KohnShamAdapter<Space>::set_poisson_rhs(const Tensor<Real, 2> &rhs)
    {
        this->poisson_rhs = rhs;
        return;
    }

    template <typename Space>
    inline void KohnShamAdapter<Space>::set_poisson_solution(const Matrix<Real> &solution)
    {
        this->poisson_solution = solution;
        this->poisson_jump_term = this->build_jump_term(solution);
        return;
    }

    template <typename Space>
    inline void KohnShamAdapter<Space>::set_kohn_sham_solution(const Matrix<Real> &solution, const List<Real> &num_elec)
    {
        this->num_elec = num_elec;
        this->kohn_sham_solution = solution;
        this->kohn_sham_jump_term = this->build_jump_term(solution);
        rho.resize(num_elec.size());
        for (Int i = 0; i < static_cast<Int>(num_elec.size()); i++)
        {
            rho[i] = num_elec[i] * this->space.get_value(kohn_sham_solution.col(i));
        }
        return;
    }

    template <typename Space>
    inline void KohnShamAdapter<Space>::set_energy(const Vector<Real> &energy)
    {
        this->energy = energy;
        return;
    }

    template <typename Space>
    inline void KohnShamAdapter<Space>::set_ext_value(const Tensor<Real, 2> &ext_value)
    {
        this->ext_value = ext_value;
        return;
    }

    template <typename Space>
    inline void KohnShamAdapter<Space>::set_har_value(const Tensor<Real, 2> &har_value)
    {
        this->har_value = har_value;
        return;
    }

    template <typename Space>
    inline void KohnShamAdapter<Space>::set_vxc_value(const Tensor<Real, 2> &vxc_value)
    {
        this->vxc_value = vxc_value;
        return;
    }

    template <typename Space>
    inline Real KohnShamAdapter<Space>::h_indicator(const Int &pid) const
    {
        const Matrix<Real, 3, 3> &inv_jacobi_matrix = this->space.get_poly(pid).inv_jacobi;
        const Matrix<Real, 3, 3> inv_jacobi_matrix_transpose = inv_jacobi_matrix.transpose();
        Tensor<Real, 2> trace(this->space.QUAD.tet.size(), this->space.INDEX_POLY_ALL.size());
        for (Int l = 0; l < static_cast<Int>(this->space.QUAD.tet.size()); l++)
        {
            for (Int it = 0; it < static_cast<Int>(this->space.INDEX_POLY_ALL.size()); it++)
            {
                trace(l, it) = ((inv_jacobi_matrix_transpose * this->space.hess_psi_value(this->space.INDEX_POLY_ALL[it])[l] * inv_jacobi_matrix).trace());
            }
        }

        Real kohn_sham_interior_residual = 0.0;
        Real kohn_sham_boundary_residual = 0.0;
        Real poisson_interior_residual = 0.0;
        Real poisson_boundary_residual = 0.0;

        for (Int i = 0; i < static_cast<Int>(num_elec.size()); i++)
        {
            for (Int l = 0; l < static_cast<Int>(this->space.QUAD.tet.size()); l++)
            {
                Real tmp = 0.0;
                for (Int it = 0; it < static_cast<Int>(this->space.INDEX_POLY_ALL.size()); it++)
                {
                    tmp -= kohn_sham_solution(this->space.num_dof_poly_all[pid][it], i) * trace(l, it);
                }
                tmp /= 2.0;
                tmp += (ext_value(pid, l) + har_value(pid, l) + vxc_value(pid, l));
                tmp -= energy(i);
                tmp *= rho[i](pid, l);
                kohn_sham_interior_residual += (this->space.get_poly(pid).volume * tmp * tmp);
            }
        }

        for (Int l = 0; l < static_cast<Int>(this->space.QUAD.tet.size()); l++)
        {
            Real tmp = 0.0;
            tmp += poisson_rhs(pid, l);
            for (Int it = 0; it < static_cast<Int>(this->space.INDEX_POLY_ALL.size()); it++)
            {
                tmp += poisson_solution(this->space.num_dof_poly_all[pid][it]) * trace(l, it);
            }
            poisson_interior_residual += (this->space.get_poly(pid).volume * tmp * tmp);
        }

        const Real h_p = this->space.get_poly_diameter(pid) / static_cast<Real>(this->space.MAX_ORDER);
        kohn_sham_interior_residual *= (h_p * h_p);
        poisson_interior_residual *= (h_p * h_p);

        for (Int i = 0; i < static_cast<Int>(num_elec.size()); i++)
        {
            for (const Int fid : this->space.get_poly(pid).face)
            {
                kohn_sham_boundary_residual += num_elec[i] * kohn_sham_jump_term(i, fid);
            }
        }
        for (const Int fid : this->space.get_poly(pid).face)
        {
            poisson_boundary_residual += poisson_jump_term(0, fid);
        }

        return std::sqrt(kohn_sham_interior_residual + kohn_sham_boundary_residual / 2.0);
    }
};

#endif
