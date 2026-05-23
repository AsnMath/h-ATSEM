#ifndef ASN_MATH_SPACE_UTILS_HPP
#define ASN_MATH_SPACE_UTILS_HPP

#include "Utils.h"

namespace Asn::Math
{
    template <typename Space>
    inline Utils<Space>::Utils(Space &space) : space(space), index_dof(_index_dof), dof_index(_dof_index) {}

    template <typename Space>
    inline void Utils<Space>::init()
    {
        _dof_index.clear();
        _index_dof.clear();
        for (Int dof = 0, index = 0; dof < space.num_dof; dof++)
        {
            if (!space.boundary_dof.contains(dof))
            {
                _index_dof[index] = dof;
                _dof_index[dof] = index;
                index++;
            }
        }
        return;
    }

    template <typename Space>
    inline Vector<Real> Utils<Space>::interpolate(const Func<Real(const Point<Space::MESH_DIM> &)> &u) const
    {
        Vector<Real> res = Vector<Real>::Zero(space.num_dof);
        #pragma omp parallel for default(none) shared(u, res)
        for (Int vid = 0; vid < space.num_vert(); vid++)
        {
            res(space.num_dof_vert[vid]) = u(space.get_vert(vid).coord);
        }
        #pragma omp parallel for default(none) shared(u, res)
        for (Int eid = 0; eid < space.num_edge(); eid++)
        {
            space.interp_edge(res, u, eid);
        }
        #pragma omp parallel for default(none) shared(u, res)
        for (Int fid = 0; fid < space.num_face(); fid++)
        {
            space.interp_face(res, u, fid);
        }
        #pragma omp parallel for default(none) shared(u, res)
        for (Int pid = 0; pid < space.num_poly(); pid++)
        {
            space.interp_poly(res, u, pid);
        }
        return res;
    }

    template <typename Space>
    inline SparseMatrix<Real> Utils<Space>::get_stiff_matrix(const MatrixMode &mode) const
    {
        if (mode == MatrixMode::Full)
        {
            return space.stiff_matrix;
        }
        else if (mode == MatrixMode::DelBoundary)
        {
            return space.map_matrix(space.stiff_matrix, dof_index);
        }
        assert(false);
    }

    template <typename Space>
    inline SparseMatrix<Real> Utils<Space>::get_mass_matrix(const MatrixMode &mode) const
    {
        if (mode == MatrixMode::Full)
        {
            return space.mass_matrix;
        }
        else if (mode == MatrixMode::DelBoundary)
        {
            return space.map_matrix(space.mass_matrix, dof_index);
        }
        assert(false);
    }
};

#endif
