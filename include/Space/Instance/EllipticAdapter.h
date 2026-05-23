#ifndef ASN_MATH_SPACE_ELLIPTIC_ADAPTER_H
#define ASN_MATH_SPACE_ELLIPTIC_ADAPTER_H

#include "../Adapter.h"

namespace Asn::Math
{
    template <typename Space>
    class EllipticAdapter;

    template <typename Space>
    class EllipticAdapter final : public Adapter<Space>
    {
    private:
        Tensor<Real, 2> rhs;
        Tensor<Real, 2> jump_term;
        Matrix<Real> solution;

    public:
        explicit EllipticAdapter(Space &space);
        EllipticAdapter(EllipticAdapter &adapter) = delete;
        virtual ~EllipticAdapter() override = default;

        void set_solution(const Matrix<Real> &solution);
        void set_rhs(const Tensor<Real, 2> &rhs);

        virtual Real h_indicator(const Int &pid) const override;

        EllipticAdapter &operator=(EllipticAdapter &adapter) = delete;
    };
};

#endif
