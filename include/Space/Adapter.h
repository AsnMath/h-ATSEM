#ifndef ASN_MATH_SPACE_ADAPTER_H
#define ASN_MATH_SPACE_ADAPTER_H

#include "../Config.h"

namespace Asn::Math
{
    template <typename Space>
    class Adapter;

    template <typename Space>
    class Adapter
    {
    protected:
        Space &space;

    public:
        explicit Adapter(Space &space);
        Adapter(Adapter &adapter) = delete;
        virtual ~Adapter() = default;

        Tensor<Real, 2> build_jump_term(const Matrix<Real> &solution) const;
        virtual Real h_indicator(const Int &pid) const = 0;
        void set_h_adaptive_flag(const String &mode, const Real &rate) const;
    };
};

#endif
