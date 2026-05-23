#ifndef ASN_MATH_SPACE_KOHN_SHAM_UTILS_HPP
#define ASN_MATH_SPACE_KOHN_SHAM_UTILS_HPP

#include <xc.h>

#include "../../../3rdparty/Eigen/QR"

#include "KohnShamUtils.h"

#include "EllipticUtils.hpp"

namespace Asn::Math
{
    template <typename Space>
    inline Real KohnShamUtils<Space>::element_effect_matrix(const Int &pid, const Index<3> &i, const Index<3> &j, const Tensor<Real, 2> &eff) const
    {
        Real res_v = 0.0;
        for (Int l = 0; l < static_cast<Int>(this->space.QUAD.tet.size()); l++)
        {
            res_v += this->space.QUAD.tet[l].second * eff(pid, l) * this->space.psi_value(i)[l] * this->space.psi_value(j)[l];
        }
        return this->space.get_poly(pid).volume * res_v;
    }

    template <typename Space>
    inline KohnShamUtils<Space>::KohnShamUtils(Space &space) : EllipticUtils<Space>(space) {}

    template <typename Space>
    inline void KohnShamUtils<Space>::init()
    {
        EllipticUtils<Space>::init();
        this->space.build_mass_matrix();
        return;
    }

    template <typename Space>
    inline Tensor<Real, 2> KohnShamUtils<Space>::get_rho(const Func<Real(const Real &)> &rho, const List<Point<Space::MESH_DIM>> &pos, const List<Real> &charge)
    const
    {
        Tensor<Real, 2> value(this->space.num_poly(), this->space.QUAD.tet.size());
        value.setZero();
        const Int num_poly = this->space.num_poly();
        #pragma omp parallel for default(none) shared(num_poly, value, rho, pos, charge)
        for (Int pid = 0; pid < num_poly; pid++)
        {
            const Matrix<Real, 3, 3> &jacobi_matrix = this->space.get_poly(pid).jacobi;
            const Point<Space::MESH_DIM> &p_0 = this->space.get_vert(this->space.get_poly(pid).vert[0]).coord;
            for (Int l = 0; l < static_cast<Int>(this->space.QUAD.tet.size()); l++)
            {
                const Point<Space::MESH_DIM> &x = p_0 + jacobi_matrix * this->space.QUAD.tet[l].first;
                for (Int i = 0; i < static_cast<Int>(pos.size()); i++)
                {
                    value(pid, l) += charge[i] * rho((x - pos[i]).template lpNorm<2>());
                }
            }
        }
        return value;
    }

    template <typename Space>
    inline Tensor<Real, 2> KohnShamUtils<Space>::get_rho(const Matrix<Real> &x, const List<Real> &F) const
    {
        assert(static_cast<Int>(F.size()) == x.cols());
        Tensor<Real, 2> value(this->space.num_poly(), this->space.QUAD.tet.size());
        value.setZero();
        for (Int i = 0; i < static_cast<Int>(F.size()); i++)
        {
            const Tensor<Real, 2> tmp = this->space.get_value(this->space.l2_normalize(x.col(i)));
            value = value + F[i] * tmp * tmp;
        }
        return value;
    }

    template <typename Space>
    inline Tensor<Real, 2> KohnShamUtils<Space>::get_ext(const List<Point<Space::MESH_DIM>> &pos, const List<Real> &F) const
    {
        const Int N = pos.size();
        Tensor<Real, 2> ext(this->space.num_poly(), this->space.QUAD.tet.size());
        for (Int pid = 0; pid < this->space.num_poly(); pid++)
        {
            const Vector<Real> &p_0 = this->space.get_vert(this->space.get_poly(pid).vert[0]).coord;
            const Matrix<Real> &jacobi_matrix = this->space.get_poly(pid).jacobi;
            for (Int l = 0; l < static_cast<Int>(this->space.QUAD.tet.size()); l++)
            {
                const Vector<Real> x = p_0 + jacobi_matrix * this->space.QUAD.tet[l].first;
                ext(pid, l) = 0.0;
                for (Int i = 0; i < N; i++)
                {
                    ext(pid, l) -= (F[i] / (x - pos[i]).template lpNorm<2>());
                }
            }
        }
        return ext;
    }

    template <typename Space>
    inline Tuple<Tensor<Real, 2>, Tensor<Real, 2>> KohnShamUtils<Space>::get_lda(const Tensor<Real, 2> &rho) const
    {
        const Int N = this->space.QUAD.tet.size();
        Tensor<Real, 2> vxc(this->space.num_poly(), this->space.QUAD.tet.size());
        Tensor<Real, 2> exc(this->space.num_poly(), this->space.QUAD.tet.size());
        #pragma omp parallel for default(none) shared(rho, vxc, exc, N)
        for (Int pid = 0; pid < this->space.num_poly(); pid++)
        {
            List<Real> r(N);
            List<Real> vx(N);
            List<Real> vc(N);
            List<Real> ex(N);
            List<Real> ec(N);
            for (Int l = 0; l < N; l++)
            {
                r[l] = rho(pid, l);
            }
            xc_func_type fx, fc;
            int tmp = xc_func_init(&fx, LIBXC_LDA_X_ID, XC_UNPOLARIZED);
            assert(tmp == 0);
            tmp = xc_func_init(&fc, LIBXC_LDA_C_ID, XC_UNPOLARIZED);
            assert(tmp == 0);
            xc_lda_exc_vxc(&fc, N, r.data(), ec.data(), vc.data());
            xc_lda_exc_vxc(&fx, N, r.data(), ex.data(), vx.data());
            xc_func_end(&fx);
            xc_func_end(&fc);
            for (Int l = 0; l < N; l++)
            {
                vxc(pid, l) = vx[l] + vc[l];
                exc(pid, l) = ex[l] + ec[l];
            }
        }
        return std::make_tuple(vxc, exc);
    }

    template <typename Space>
    inline Tensor<Real, 2> KohnShamUtils<Space>::get_lda_vxc(const Tensor<Real, 2> &rho) const
    {
        const Int N = this->space.num_poly() * this->space.QUAD.tet.size();
        List<Real> r(N);
        List<Real> vx(N);
        List<Real> vc(N);
        for (Int pid = 0, i = 0; pid < this->space.num_poly(); pid++)
        {
            for (Int l = 0; l < static_cast<Int>(this->space.QUAD.tet.size()); l++, i++)
            {
                r[i] = rho(pid, l);
            }
        }
        xc_func_type fx, fc;
        int tmp = xc_func_init(&fx, LIBXC_LDA_X_ID, XC_UNPOLARIZED);
        assert(tmp == 0);
        xc_lda_vxc(&fx, N, r.data(), vx.data());
        xc_func_end(&fx);
        tmp = xc_func_init(&fc, LIBXC_LDA_C_ID, XC_UNPOLARIZED);
        assert(tmp == 0);
        xc_lda_vxc(&fc, N, r.data(), vc.data());
        xc_func_end(&fc);
        Tensor<Real, 2> vxc(this->space.num_poly(), this->space.QUAD.tet.size());
        for (Int pid = 0, i = 0; pid < this->space.num_poly(); pid++)
        {
            for (Int l = 0; l < static_cast<Int>(this->space.QUAD.tet.size()); l++, i++)
            {
                vxc(pid, l) = vx[i] + vc[i];
            }
        }
        return vxc;
    }

    template <typename Space>
    inline Tensor<Real, 2> KohnShamUtils<Space>::get_lda_exc(const Tensor<Real, 2> &rho) const
    {
        const Int N = this->space.num_poly() * this->space.QUAD.tet.size();
        List<Real> r(N);
        List<Real> ex(N);
        List<Real> ec(N);
        for (Int pid = 0, i = 0; pid < this->space.num_poly(); pid++)
        {
            for (Int l = 0; l < static_cast<Int>(this->space.QUAD.tet.size()); l++, i++)
            {
                r[i] = rho(pid, l);
            }
        }
        xc_func_type fx, fc;
        int tmp = xc_func_init(&fx, LIBXC_LDA_X_ID, XC_UNPOLARIZED);
        assert(tmp == 0);
        xc_lda_exc(&fx, N, r.data(), ex.data());
        xc_func_end(&fx);
        tmp = xc_func_init(&fc, LIBXC_LDA_C_ID, XC_UNPOLARIZED);
        assert(tmp == 0);
        xc_lda_exc(&fc, N, r.data(), ec.data());
        xc_func_end(&fc);
        Tensor<Real, 2> exc(this->space.num_poly(), this->space.QUAD.tet.size());
        for (Int pid = 0, i = 0; pid < this->space.num_poly(); pid++)
        {
            for (Int l = 0; l < static_cast<Int>(this->space.QUAD.tet.size()); l++, i++)
            {
                exc(pid, l) = ex[i] + ec[i];
            }
        }
        return exc;
    }

    template <typename Space>
    inline Tuple<Tensor<Real, 2>, Tensor<Real, 2>, Tensor<Real, 2>> KohnShamUtils<Space>::get_lsd(const Tensor<Real, 2> &rho_1, const Tensor<Real, 2> &rho_2) const
    {
        const Int N = this->space.QUAD.tet.size();
        Tensor<Real, 2> vxc_1(this->space.num_poly(), this->space.QUAD.tet.size());
        Tensor<Real, 2> vxc_2(this->space.num_poly(), this->space.QUAD.tet.size());
        Tensor<Real, 2> exc(this->space.num_poly(), this->space.QUAD.tet.size());
        #pragma omp parallel for default(none) shared(rho_1, rho_2, vxc_1, vxc_2, exc, N)
        for (Int pid = 0; pid < this->space.num_poly(); pid++)
        {
            List<Real> r(2 * N);
            List<Real> vx(2 * N);
            List<Real> vc(2 * N);
            List<Real> ex(N);
            List<Real> ec(N);
            for (Int l = 0; l < static_cast<Int>(this->space.QUAD.tet.size()); l++)
            {
                r[2 * l] = rho_1(pid, l);
                r[2 * l + 1] = rho_2(pid, l);
            }
            xc_func_type fx, fc;
            int tmp = xc_func_init(&fx, LIBXC_LDA_X_ID, XC_POLARIZED);
            assert(tmp == 0);
            tmp = xc_func_init(&fc, LIBXC_LDA_C_ID, XC_POLARIZED);
            assert(tmp == 0);
            xc_lda_exc_vxc(&fc, N, r.data(), ec.data(), vc.data());
            xc_lda_exc_vxc(&fx, N, r.data(), ex.data(), vx.data());
            xc_func_end(&fx);
            xc_func_end(&fc);
            for (Int l = 0; l < static_cast<Int>(this->space.QUAD.tet.size()); l++)
            {
                vxc_1(pid, l) = vx[2 * l] + vc[2 * l];
                vxc_2(pid, l) = vx[2 * l + 1] + vc[2 * l + 1];
                exc(pid, l) = ex[l] + ec[l];
            }
        }
        return std::make_tuple(vxc_1, vxc_2, exc);
    }

    template <typename Space>
    inline Pair<Tensor<Real, 2>, Tensor<Real, 2>> KohnShamUtils<Space>::get_lsd_vxc(const Tensor<Real, 2> &rho_1, const Tensor<Real, 2> &rho_2) const
    {
        const Int N = this->space.num_poly() * this->space.QUAD.tet.size();
        List<Real> r(2 * N);
        List<Real> vx(2 * N);
        List<Real> vc(2 * N);
        for (Int pid = 0, i = 0; pid < this->space.num_poly(); pid++)
        {
            for (Int l = 0; l < static_cast<Int>(this->space.QUAD.tet.size()); l++, i += 2)
            {
                r[i] = rho_1(pid, l);
                r[i + 1] = rho_2(pid, l);
            }
        }
        xc_func_type fx, fc;
        int tmp = xc_func_init(&fx, LIBXC_LDA_X_ID, XC_POLARIZED);
        assert(tmp == 0);
        xc_lda_vxc(&fx, N, r.data(), vx.data());
        xc_func_end(&fx);
        tmp = xc_func_init(&fc, LIBXC_LDA_C_ID, XC_POLARIZED);
        assert(tmp == 0);
        xc_lda_vxc(&fc, N, r.data(), vc.data());
        xc_func_end(&fc);
        Tensor<Real, 2> vxc_1(this->space.num_poly(), this->space.QUAD.tet.size());
        Tensor<Real, 2> vxc_2(this->space.num_poly(), this->space.QUAD.tet.size());
        for (Int pid = 0, i = 0; pid < this->space.num_poly(); pid++)
        {
            for (Int l = 0; l < static_cast<Int>(this->space.QUAD.tet.size()); l++, i += 2)
            {
                vxc_1(pid, l) = vx[i] + vc[i];
                vxc_2(pid, l) = vx[i + 1] + vc[i + 1];
            }
        }
        return std::make_pair(vxc_1, vxc_2);
    }

    template <typename Space>
    inline Tensor<Real, 2> KohnShamUtils<Space>::get_lsd_exc(const Tensor<Real, 2> &rho_1, const Tensor<Real, 2> &rho_2) const
    {
        const Int N = this->space.num_poly() * this->space.QUAD.tet.size();
        List<Real> r(2 * N);
        List<Real> ex(N);
        List<Real> ec(N);
        for (Int pid = 0, i = 0; pid < this->space.num_poly(); pid++)
        {
            for (Int l = 0; l < static_cast<Int>(this->space.QUAD.tet.size()); l++, i += 2)
            {
                r[i] = rho_1(pid, l);
                r[i + 1] = rho_2(pid, l);
            }
        }
        xc_func_type fx, fc;
        int tmp = xc_func_init(&fx, LIBXC_LDA_X_ID, XC_POLARIZED);
        assert(tmp == 0);
        xc_lda_exc(&fx, N, r.data(), ex.data());
        xc_func_end(&fx);
        tmp = xc_func_init(&fc, LIBXC_LDA_C_ID, XC_POLARIZED);
        assert(tmp == 0);
        xc_lda_exc(&fc, N, r.data(), ec.data());
        xc_func_end(&fc);
        Tensor<Real, 2> exc(this->space.num_poly(), this->space.QUAD.tet.size());
        for (Int pid = 0, i = 0; pid < this->space.num_poly(); pid++)
        {
            for (Int l = 0; l < static_cast<Int>(this->space.QUAD.tet.size()); l++, i++)
            {
                exc(pid, l) = ex[i] + ec[i];
            }
        }
        return exc;
    }

    template <typename Space>
    inline SparseMatrix<Real> KohnShamUtils<Space>::get_effect_matrix(const Tensor<Real, 2> &eff) const
    {
        List<List<Eigen::Triplet<Real>>> list(this->space.num_poly());
        #pragma omp parallel for default(none) shared(list, eff)
        for (Int pid = 0; pid < this->space.num_poly(); pid++)
        {
            for (Int it = 0; it < static_cast<Int>(this->space.INDEX_POLY_ALL.size()); it++)
            {
                const Int dof_1 = this->space.num_dof_poly_all[pid][it];
                if (!this->dof_index.contains(dof_1))
                {
                    continue;
                }
                const Int index_1 = this->dof_index.at(dof_1);
                list[pid].emplace_back(index_1, index_1, element_effect_matrix(pid, this->space.INDEX_POLY_ALL[it], this->space.INDEX_POLY_ALL[it], eff));
                for (Int jt = it + 1; jt < static_cast<Int>(this->space.INDEX_POLY_ALL.size()); jt++)
                {
                    const Int dof_2 = this->space.num_dof_poly_all[pid][jt];
                    if (!this->dof_index.contains(dof_2))
                    {
                        continue;
                    }
                    if (const Int index_2 = this->dof_index.at(dof_2); index_1 > index_2)
                    {
                        list[pid].emplace_back(index_1, index_2, element_effect_matrix(pid, this->space.INDEX_POLY_ALL[it], this->space.INDEX_POLY_ALL[jt], eff));
                    }
                    else
                    {
                        list[pid].emplace_back(index_2, index_1, element_effect_matrix(pid, this->space.INDEX_POLY_ALL[it], this->space.INDEX_POLY_ALL[jt], eff));
                    }
                }
            }
        }
        List<Eigen::Triplet<Real>> list_all;
        Int size = 0;
        for (Int pid = 0; pid < this->space.num_poly(); pid++)
        {
            size += list[pid].size();
        }
        list_all.reserve(size);
        for (Int pid = 0; pid < this->space.num_poly(); pid++)
        {
            list_all.insert(list_all.end(), list[pid].begin(), list[pid].end());
        }
        SparseMatrix<Real> A(this->index_dof.size(), this->index_dof.size());
        A.setFromTriplets(list_all.begin(), list_all.end());
        A.prune(static_cast<Real>(0.0));
        const SparseMatrix<Real> effect_matrix = A.selfadjointView<Eigen::Lower>();
        return effect_matrix;
    }

    inline Real delta(const Int &i, const Int &j)
    {
        return i == j ? 1.0 : 0.0;
    }

    template <typename Space>
    inline Tuple<Real, Vector<Real>, Matrix<Real>, Tensor<Real, 3>> KohnShamUtils<Space>::get_multipole_expansion(const Tensor<Real, 2> &rho) const
    {
        const Int N = this->space.num_poly();
        Real q = 0.0;
        Vector<Real> p = Vector<Real>::Zero(3);
        Matrix<Real> Q = Matrix<Real>::Zero(3, 3);
        Tensor<Real, 3> O(3, 3, 3);
        O.setZero();
        List<Real> list_q(N);
        List<Vector<Real>> list_p(N);
        List<Matrix<Real>> list_Q(N);
        List<Tensor<Real, 3>> list_O(N);
        const Matrix<Real> I = Matrix<Real>::Identity(3, 3);
        #pragma omp parallel for default(none) shared(list_q, list_p, list_Q, list_O, I, rho, N)
        for (Int pid = 0; pid < N; pid++)
        {
            Real tmp_q = 0.0;
            Vector<Real> tmp_p = Vector<Real>::Zero(3);
            Matrix<Real> tmp_Q = Matrix<Real>::Zero(3, 3);
            Tensor<Real, 3> tmp_O(3, 3, 3);
            tmp_O.setZero();
            const Vector<Real> &p_0 = this->space.get_vert(this->space.get_poly(pid).vert[0]).coord;
            const Matrix<Real> &jacobi_matrix = this->space.get_poly(pid).jacobi;
            for (Int l = 0; l < static_cast<Int>(this->space.QUAD.tet.size()); l++)
            {
                const Vector<Real> x = p_0 + jacobi_matrix * this->space.QUAD.tet[l].first;
                const Real squared_norm = x.squaredNorm();
                tmp_q += this->space.QUAD.tet[l].second * rho(pid, l);
                tmp_p += this->space.QUAD.tet[l].second * rho(pid, l) * x;
                tmp_Q += this->space.QUAD.tet[l].second * rho(pid, l) * (3.0 * x * x.transpose() - x.squaredNorm() * I);
                for (Int i = 0; i < 3; i++)
                {
                    for (Int j = 0; j < 3; j++)
                    {
                        for (Int k = 0; k < 3; k++)
                        {
                            tmp_O(i, j, k) += this->space.QUAD.tet[l].second * rho(pid, l) * (15.0 * x(i) * x(j) * x(k) - 3.0 * squared_norm * (x(i) * delta(j, k) + x(j) * delta(i, k) + x(k) * delta(i, j)));
                        }
                    }
                }
            }
            const Real volume = this->space.get_poly(pid).volume;
            list_q[pid] = tmp_q * volume;
            list_p[pid] = tmp_p * volume;
            list_Q[pid] = tmp_Q * volume;
            list_O[pid] = tmp_O * volume;
        }
        for (Int pid = 0; pid < N; pid++)
        {
            q += list_q[pid];
            p += list_p[pid];
            Q += list_Q[pid];
            O += list_O[pid];
        }
        return std::make_tuple(q, p, Q, O);
    }

    template <typename Space>
    inline Real KohnShamUtils<Space>::get_energy_ion_ion(const List<Point<Space::MESH_DIM>> &pos, const List<Real> &charge)
    {
        Real res = 0.0;
        const Int N = pos.size();
        for (Int i = 0; i < N; i++)
        {
            for (Int j = i + 1; j < N; j++)
            {
                res += (charge[i] * charge[j] / (pos[i] - pos[j]).template lpNorm<2>());
            }
        }
        return res;
    }

    template <typename Space>
    inline Vector<Real> KohnShamUtils<Space>::get_energy(const Vector<Real> &eps, const List<Real> &F, const Matrix<Real> &solution,
                                                         const Tensor<Real, 2> &har, const Tensor<Real, 2> &vxc, const Tensor<Real, 2> &exc) const
    {
        Vector<Real> energy(F.size());
        const Tensor<Real, 2> value = exc - vxc - har / 2.0;
        List<Real> list(this->space.num_poly());
        for (Int i = 0; i < static_cast<Int>(F.size()); i++)
        {
            Tensor<Real, 2> rho = get_rho(solution.col(i), {F[i]});
            #pragma omp parallel for default(none) shared(value, list, rho)
            for (Int pid = 0; pid < this->space.num_poly(); pid++)
            {
                list[pid] = 0.0;
                for (Int l = 0; l < static_cast<Int>(this->space.QUAD.tet.size()); l++)
                {
                    list[pid] += this->space.QUAD.tet[l].second * value(pid, l) * rho(pid, l);
                }
                list[pid] *= this->space.get_poly(pid).volume;
            }
            energy(i) = F[i] * eps(i);
            for (Int pid = 0; pid < this->space.num_poly(); pid++)
            {
                energy(i) += list[pid];
            }
        }
        return energy;
    }

    template <typename Space>
    inline Tuple<Vector<Real>, Vector<Real>> KohnShamUtils<Space>::get_energy(const Vector<Real> &eps_up, const Vector<Real> &eps_down,
                                                                              const List<Real> &F_up, const List<Real> &F_down,
                                                                              const Matrix<Real> &solution_up, const Matrix<Real> &solution_down, const Tensor<Real, 2> &har,
                                                                              const Tensor<Real, 2> &vxc_up, const Tensor<Real, 2> &vxc_down, const Tensor<Real, 2> &exc) const
    {
        Vector<Real> energy_up(F_up.size());
        Vector<Real> energy_down(F_down.size());
        const Tensor<Real, 2> value_up = exc - har / 2.0 - vxc_up;
        const Tensor<Real, 2> value_down = exc - har / 2.0 - vxc_down;
        List<Real> list(this->space.num_poly());
        for (Int i = 0; i < static_cast<Int>(F_up.size()); i++)
        {
            const Tensor<Real, 2> rho = get_rho(solution_up.col(i), {F_up[i]});
            #pragma omp parallel for default(none) shared(list, solution_up, F_up, i, value_up, rho)
            for (Int pid = 0; pid < this->space.num_poly(); pid++)
            {
                list[pid] = 0.0;
                for (Int l = 0; l < static_cast<Int>(this->space.QUAD.tet.size()); l++)
                {
                    list[pid] += (this->space.QUAD.tet[l].second * value_up(pid, l) * rho(pid, l));
                }
                list[pid] *= this->space.get_poly(pid).volume;
            }
            energy_up(i) = F_up[i] * eps_up(i);
            for (Int pid = 0; pid < this->space.num_poly(); pid++)
            {
                energy_up(i) += list[pid];
            }
        }
        for (Int i = 0; i < static_cast<Int>(F_down.size()); i++)
        {
            const Tensor<Real, 2> rho = get_rho(solution_down.col(i), {F_down[i]});
            #pragma omp parallel for default(none) shared(list, solution_down, F_down, i, value_down, rho)
            for (Int pid = 0; pid < this->space.num_poly(); pid++)
            {
                list[pid] = 0.0;
                for (Int l = 0; l < static_cast<Int>(this->space.QUAD.tet.size()); l++)
                {
                    list[pid] += (this->space.QUAD.tet[l].second * value_down(pid, l) * rho(pid, l));
                }
                list[pid] *= this->space.get_poly(pid).volume;
            }
            energy_down(i) = F_down[i] * eps_down(i);
            for (Int pid = 0; pid < this->space.num_poly(); pid++)
            {
                energy_down(i) += list[pid];
            }
        }
        return std::make_tuple(energy_up, energy_down);
    }

    template <typename Space>
    inline Tensor<Real, 2> KohnShamUtils<Space>::pulay_mixing(const List<Tensor<Real, 2>> &rho, const List<Tensor<Real, 2>> &res, const Real &factor) const
    {
        const Int N = rho.size();
        Matrix<Real> A(N + 1, N + 1);
        A(N, N) = 0.0;
        for (Int i = 0; i < N; i++)
        {
            A(i, i) = this->space.inner_product(res[i], res[i]);
            for (Int j = i + 1; j < N; j++)
            {
                A(i, j) = A(j, i) = this->space.inner_product(res[i], res[j]);
            }
            A(i, N) = 1.0;
            A(N, i) = 1.0;
        }
        Vector<Real> rhs = Vector<Real>::Zero(N + 1);
        rhs(N) = 1.0;
        Vector<Real> alpha = A.completeOrthogonalDecomposition().solve(rhs);
        Tensor<Real, 2> rho_new(rho[0].dimensions());
        rho_new.setZero();
        for (Int i = 0; i < N; i++)
        {
            rho_new += alpha(i) * (rho[i] + factor * res[i]);
        }
        return rho_new;
    }
};

#endif
