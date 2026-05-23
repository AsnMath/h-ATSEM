#ifndef ASN_MATH_SPACE_SPLINE_HPP
#define ASN_MATH_SPACE_SPLINE_HPP

#include "Spline.h"

#include "Space.hpp"

namespace Asn::Math
{
    template <Real TOL>
    inline Bool Spline<TOL>::inside_poly(const Point<BaseMesh::MESH_DIM> &p, const Array<Point<BaseMesh::MESH_DIM>, 4> &t)
    {
        for (Int i = 0; i < BaseMesh::MESH_DIM; i++)
        {
            const Real m = std::min(std::min(t[0](i), t[1](i)), std::min(t[2](i), t[3](i))) - TOL;
            const Real M = std::max(std::max(t[0](i), t[1](i)), std::max(t[2](i), t[3](i))) + TOL;
            if (p(i) < m || p(i) > M)
            {
                return false;
            }
        }
        Matrix<Real, 3, 3> jacobi_matrix;
        jacobi_matrix(0, 0) = t[1](0) - t[0](0);
        jacobi_matrix(0, 1) = t[2](0) - t[0](0);
        jacobi_matrix(0, 2) = t[3](0) - t[0](0);
        jacobi_matrix(1, 0) = t[1](1) - t[0](1);
        jacobi_matrix(1, 1) = t[2](1) - t[0](1);
        jacobi_matrix(1, 2) = t[3](1) - t[0](1);
        jacobi_matrix(2, 0) = t[1](2) - t[0](2);
        jacobi_matrix(2, 1) = t[2](2) - t[0](2);
        jacobi_matrix(2, 2) = t[3](2) - t[0](2);
        const Point<BaseMesh::MESH_DIM> rp = jacobi_matrix.fullPivLu().solve(p - t[0]);
        if (rp(0) >= -TOL && rp(1) >= -TOL && rp(2) >= -TOL && rp.sum() <= 1.0 + TOL)
        {
            return true;
        }
        return false;
    }

    template <Real TOL>
    inline Bool Spline<TOL>::inside_poly(const List<Point<BaseMesh::MESH_DIM>> &p, const Array<Point<BaseMesh::MESH_DIM>, 4> &t)
    {
        for (Int i = 0; i < static_cast<Int>(p.size()); i++)
        {
            if (inside_poly(p[i], t) == false)
            {
                return false;
            }
        }
        return true;
    }

    template <Real TOL>
    inline void Spline<TOL>::init(const Space<TOL> &space)
    {
        this->BaseMesh::operator=(space);
        max_order = space.MAX_ORDER;
        index_face = space.INDEX_FACE;
        index_poly = space.INDEX_POLY;
        index_poly_all = space.INDEX_POLY_ALL;
        num_dof_vert = space.num_dof_vert;
        num_dof_edge = space.num_dof_edge;
        num_dof_face = space.num_dof_face;
        num_dof_poly = space.num_dof_poly;
        num_dof_poly_all = space.num_dof_poly_all;

        return;
    }

    template <Real TOL>
    inline void Spline<TOL>::init(const Matrix<Real> &solution)
    {
        this->solution = solution;
        return;
    }

    template <Real TOL>
    inline Real Spline<TOL>::get_value(const Point<BaseMesh::MESH_DIM> &p, const Int &index, ConstNodePtr ptr) const
    {
        ptr = find_leaf({p}, ptr);
        const Int pid = this->find_poly(ptr->data);
        const Matrix<Real, 3, 3> &inv_jacobi_matrix = this->get_poly(pid).inv_jacobi;
        const Point<BaseMesh::MESH_DIM> &p_0 = vert[this->get_poly(pid).vert[0]].coord;
        const Point<BaseMesh::MESH_DIM> rp = inv_jacobi_matrix * (p - p_0);
        Real value = 0.0;
        for (Int i = 0; i < static_cast<Int>(index_poly_all.size()); i++)
        {
            value += (solution(num_dof_poly_all[pid][i], index) * Space<TOL>::template psi<Real>(index_poly_all[i], rp));
        }
        return value;
    }

    template <Real TOL>
    inline Vector<Real> Spline<TOL>::get_value(const Point<BaseMesh::MESH_DIM> &p, ConstNodePtr ptr) const
    {
        ptr = find_leaf({p}, ptr);
        const Int pid = this->find_poly(ptr->data);
        const Matrix<Real, 3, 3> &inv_jacobi_matrix = this->get_poly(pid).inv_jacobi;
        const Point<BaseMesh::MESH_DIM> &p_0 = vert[this->get_poly(pid).vert[0]].coord;
        const Point<BaseMesh::MESH_DIM> rp = inv_jacobi_matrix * (p - p_0);
        Vector<Real> value = Vector<Real>::Zero(solution.cols());
        for (Int i = 0; i < static_cast<Int>(index_poly_all.size()); i++)
        {
            value += (solution.row(num_dof_poly_all[pid][i]) * Space<TOL>::template psi<Real>(index_poly_all[i], rp));
        }
        return value;
    }

    template <Real TOL>
    inline ConstNodePtr Spline<TOL>::find_node(const List<Point<BaseMesh::MESH_DIM>> &p, ConstNodePtr ptr) const
    {
        if (ptr == nullptr)
        {
            for (Int i = 0; i < static_cast<Int>(tree.root.size()); i++)
            {
                for (Int j = 0; j < static_cast<Int>(p.size()); j++)
                {
                    if (inside_poly(p, tree.root[i]->data))
                    {
                        ptr = tree.root[i];
                        break;
                    }
                }
            }
        }
        assert(ptr != nullptr);
        while (ptr->child.size() != 0)
        {
            for (Int i = 0; i <= static_cast<Int>(ptr->child.size()); i++)
            {
                if (i >= static_cast<Int>(ptr->child.size()))
                {
                    return ptr;
                }
                if (inside_poly(p, ptr->child[i]->data))
                {
                    ptr = ptr->child[i];
                    break;
                }
            }
        }
        return ptr;
    }

    template <Real TOL>
    inline ConstNodePtr Spline<TOL>::find_leaf(const List<Point<BaseMesh::MESH_DIM>> &p, ConstNodePtr ptr) const
    {
        ptr = find_node(p, ptr);
        if (ptr->child.size() == 0)
        {
            return ptr;
        }
        assert(false);
    }

    template <Real TOL>
    inline Matrix<Real> Spline<TOL>::interp(const Matrix<Real> &solution, const Space<TOL> &space)
    {
        this->solution = solution;
        const Int M = solution.cols();
        Matrix<Real> res = Matrix<Real>::Zero(space.num_dof, M);
        List<ConstNodePtr> ptr(space.num_poly(), nullptr);
        List<Int> id(space.num_poly(), NPOS);
        #pragma omp parallel for default(none) shared(space, ptr, id)
        for (Int pid = 0; pid < space.num_poly(); pid++)
        {
            const Array<Int, VERTS_PER_POLY> vid = space.get_poly(pid).vert;
            const Array<Point<BaseMesh::MESH_DIM>, VERTS_PER_POLY> list = {
                        space.get_vert(vid[0]).coord, space.get_vert(vid[1]).coord,
                        space.get_vert(vid[2]).coord, space.get_vert(vid[3]).coord
                    };
            id[pid] = this->find_poly(list);
            if (id[pid] != NPOS)
            {
                ptr[pid] = this->get_poly(id[pid]).ptr;
            }
            else
            {
                ptr[pid] = find_node({list[0], list[1], list[2], list[3]});
            }
        }
        #pragma omp parallel for default(none) shared(space, res, ptr, solution, M)
        for (Int vid = 0; vid < space.num_vert(); vid++)
        {
            if (const Int svid = this->find_vert(space.get_vert(vid).coord); svid != NPOS)
            {
                res.row(space.num_dof_vert[vid]) = solution.row(num_dof_vert[svid]);
            }
            else
            {
                for (Int index = 0; index < M; index++)
                {
                    Int pid = *space.get_vert(vid).poly.begin();
                    res(space.num_dof_vert[vid], index) = get_value(space.get_vert(vid).coord, index, ptr[pid]);
                }
            }
        }
        #pragma omp parallel for default(none) shared(space, res, ptr, solution, M)
        for (Int eid = 0; eid < space.num_edge(); eid++)
        {
            const Array<Int, VERTS_PER_EDGE> vid = space.get_edge(eid).vert;
            if (const Int seid = this->find_edge({space.get_vert(vid[0]).coord, space.get_vert(vid[1]).coord}); seid != NPOS)
            {
                for (Int i = 2; i <= max_order; i++)
                {
                    res.row(space.num_dof_edge[eid][i]) = solution.row(num_dof_edge[seid][i]);
                }
            }
            else
            {
                for (Int index = 0; index < M; index++)
                {
                    Eigen::Ref<Vector<Real>> col = res.col(index);
                    Int pid = *space.get_edge(eid).poly.begin();
                    Func<Real(const Point<BaseMesh::MESH_DIM> &)>
                            u = [&](const Point<BaseMesh::MESH_DIM> &p) { return get_value(p, index, ptr[pid]); };
                    space.interp_edge(col, u, eid);
                }
            }
        }
        #pragma omp parallel for default(none) shared(space, res, ptr, solution, M)
        for (Int fid = 0; fid < space.num_face(); fid++)
        {
            const Array<Int, VERTS_PER_FACE> vid = space.get_face(fid).vert;
            const Array<Point<BaseMesh::MESH_DIM>, VERTS_PER_FACE> list = {
                        space.get_vert(vid[0]).coord, space.get_vert(vid[1]).coord, space.get_vert(vid[2]).coord
                    };
            if (const Int sfid = this->find_face(list); sfid != NPOS)
            {
                for (auto &it : index_face)
                {
                    res.row(space.num_dof_face[fid](it)) = solution.row(num_dof_face[sfid](it));
                }
            }
            else
            {
                for (Int index = 0; index < M; index++)
                {
                    Eigen::Ref<Vector<Real>> col = res.col(index);
                    Int pid = *space.get_face(fid).poly.begin();
                    Func<Real(const Point<BaseMesh::MESH_DIM> &)> u = [&](const Point<BaseMesh::MESH_DIM> &p) { return get_value(p, index, ptr[pid]); };
                    space.interp_face(col, u, fid);
                }
            }
        }
        #pragma omp parallel for default(none) shared(space, res, ptr, id, solution,M)
        for (Int pid = 0; pid < space.num_poly(); pid++)
        {
            if (id[pid] != NPOS)
            {
                for (auto &it : index_poly)
                {
                    res.row(space.num_dof_poly[pid](it)) = solution.row(num_dof_poly[id[pid]](it));
                }
            }
            else
            {
                for (Int index = 0; index < M; index++)
                {
                    Eigen::Ref<Vector<Real>> col = res.col(index);
                    Func<Real(const Point<BaseMesh::MESH_DIM> &)>
                            u = [&](const Point<BaseMesh::MESH_DIM> &p) { return get_value(p, index, ptr[pid]); };
                    space.interp_poly(col, u, pid);
                }
            }
        }
        return res;
    }
};

#endif
