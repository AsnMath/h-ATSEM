#ifndef ASN_MATH_SPACE_ADAPTER_HPP
#define ASN_MATH_SPACE_ADAPTER_HPP

#include <queue>

#include "Adapter.h"

namespace Asn::Math
{
    template <typename Space>
    inline Adapter<Space>::Adapter(Space &space) : space(space) {}

    template <typename Space>
    inline Tensor<Real, 2> Adapter<Space>::build_jump_term(const Matrix<Real> &solution) const
    {
        Tensor<Real, 2> jump_term(solution.cols(), space.num_face());
        for (Int col = 0; col < solution.cols(); col++)
        {
            #pragma omp parallel for default(none) shared(col, jump_term, solution)
            for (Int fid = 0; fid < space.num_face(); fid++)
            {
                const Real h_f = space.get_face_diameter(fid) / static_cast<Real>(space.MAX_ORDER);
                const List<Int> pid = convert<HashSet<Int>, List<Int>>(space.get_face(fid).poly);
                jump_term(col, fid) = 0.0;
                if (pid.size() == 2)
                {
                    const Array<Int, Space::VERTS_PER_FACE> vid = space.get_face(fid).vert;

                    const Point<3> &v_0 = space.get_vert(vid[0]).coord;
                    const Point<3> v_1 = space.get_vert(vid[1]).coord - v_0;
                    const Point<3> v_2 = space.get_vert(vid[2]).coord - v_0;

                    const Real area = std::abs(v_1.cross(v_2).template lpNorm<2>() / 2.0);
                    const Point<3> norm = v_1.cross(v_2) / (2.0 * area);

                    const Matrix<Real, 3, 3> inv_jacobi_transpose_0 = space.get_poly(pid[0]).inv_jacobi.transpose();
                    const Matrix<Real, 3, 3> inv_jacobi_transpose_1 = space.get_poly(pid[1]).inv_jacobi.transpose();

                    Int face_index_0 = -1;
                    Int face_index_1 = -1;

                    const Array<Int, Space::VERTS_PER_POLY> &vid_0 = space.get_poly(pid[0]).vert;
                    const Array<Int, Space::VERTS_PER_POLY> &vid_1 = space.get_poly(pid[1]).vert;

                    for (Int i = 0; i < space.VERTS_PER_POLY; i++)
                    {
                        if (vid[0] != vid_0[i] && vid[1] != vid_0[i] && vid[2] != vid_0[i])
                        {
                            assert(face_index_0 == -1);
                            face_index_0 = i;
                        }
                        if (vid[0] != vid_1[i] && vid[1] != vid_1[i] && vid[2] != vid_1[i])
                        {
                            assert(face_index_1 == -1);
                            face_index_1 = i;
                        }
                    }
                    assert(face_index_0 != -1);
                    assert(face_index_1 != -1);
                    for (Int l = 0; l < static_cast<Int>(space.QUAD.tri.size()); l++)
                    {
                        Point<3> tmp = Point<3>::Zero();
                        for (const Index<3> &index : space.INDEX_POLY_ALL)
                        {
                            tmp += (solution(space.get_num_dof_poly(pid[0], index), col) * (inv_jacobi_transpose_0 * space.grad_psi_value_face(index)[face_index_0][l])).dot(norm) * norm;
                            tmp -= (solution(space.get_num_dof_poly(pid[1], index), col) * (inv_jacobi_transpose_1 * space.grad_psi_value_face(index)[face_index_1][l])).dot(norm) * norm;
                        }
                        jump_term(col, fid) += area * space.QUAD.tri[l].second * tmp.dot(tmp);
                    }
                    jump_term(col, fid) *= h_f;
                }
            }
        }
        return jump_term;
    }

    template <typename Space>
    inline void Adapter<Space>::set_h_adaptive_flag(const String &mode, const Real &rate) const
    {
        if (mode == "all")
        {
            #pragma omp parallel for default(none)
            for (Int pid = 0; pid < space.num_poly(); pid++)
            {
                space.set_refine_flag(pid);
            }
            return;
        }
        if (mode == "none")
        {
            return;
        }

        List<Real> residual(space.num_poly());
        #pragma omp parallel for default(none) shared(residual)
        for (Int pid = 0; pid < space.num_poly(); pid++)
        {
            residual[pid] = this->h_indicator(pid);
        }

        Real refine_threshold = 0.0;

        if (mode == "rate_max")
        {
            Real max_error = 0.0;
            for (const Real &item : residual)
            {
                assert(item >= 0.0);
                if (item > max_error)
                {
                    max_error = item;
                }
            }
            refine_threshold = max_error * rate;
        }
        else if (mode == "rate_acc")
        {
            Real s = 0.0;
            for (const Real &item : residual)
            {
                assert(item >= 0.0);
                s += item;
            }
            std::priority_queue<Real> heap(residual.begin(), residual.end());
            s *= rate;
            Real acc = 0.0;
            while (!heap.empty())
            {
                acc += heap.top();
                if (acc >= s)
                {
                    refine_threshold = heap.top();
                    break;
                }
                heap.pop();
            }
        }
        else
        {
            assert(false);
        }

        #pragma omp parallel for default(none) shared(residual, refine_threshold)
        for (Int pid = 0; pid < space.num_poly(); pid++)
        {
            if (residual[pid] >= refine_threshold)
            {
                space.set_refine_flag(pid);
            }
        }
        return;
    }
};

#endif
