#ifndef ASN_MATH_SPACE_UTILS_H
#define ASN_MATH_SPACE_UTILS_H

#include "../Config.h"

namespace Asn::Math
{
    enum class MatrixMode
    {
        Full,
        DelBoundary,
    };

    template <typename Space>
    class Utils;

    template <typename Space>
    class Utils
    {
    protected:
        Space &space;

        HashMap<Int, Int> _index_dof;
        HashMap<Int, Int> _dof_index;

    public:
        const HashMap<Int, Int> &index_dof;
        const HashMap<Int, Int> &dof_index;

    public:
        explicit Utils(Space &space);
        Utils(const Utils &abstract_utils) = delete;
        virtual ~Utils() = default;

        virtual void init();

        Vector<Real> interpolate(const Func<Real(const Point<Space::MESH_DIM> &)> &u) const;

        SparseMatrix<Real> get_stiff_matrix(const MatrixMode &mode) const;
        SparseMatrix<Real> get_mass_matrix(const MatrixMode &mode) const;

        Utils &operator=(const Utils &abstract_utils) = delete;
    };
};

#endif
