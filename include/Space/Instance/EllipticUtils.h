#ifndef ASN_MATH_SPACE_ELLIPTIC_UTILS_H
#define ASN_MATH_SPACE_ELLIPTIC_UTILS_H

#include "../Utils.h"

namespace Asn::Math
{
    template <typename Space>
    class EllipticUtils;

    template <typename Space>
    class EllipticUtils : public Utils<Space>
    {
    public:
        explicit EllipticUtils(Space &space);
        EllipticUtils(const EllipticUtils &abstract_utils) = delete;
        virtual ~EllipticUtils() override = default;

        virtual void init() override;

        Vector<Real> get_boundary(const Func<Real(const Point<Space::MESH_DIM> &)> &u) const;
        Vector<Real> get_rhs(const Tensor<Real, 2> &f) const;
        void apply_boundary(SparseMatrix<Real> &A, Vector<Real> &u) const;

        EllipticUtils &operator=(const EllipticUtils &abstract_utils) = delete;
    };
};

#endif
