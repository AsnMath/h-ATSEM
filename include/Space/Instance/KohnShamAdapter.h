#ifndef ASN_MATH_SPACE_KOHN_SHAM_ADAPTER_H
#define ASN_MATH_SPACE_KOHN_SHAM_ADAPTER_H

#include "../Adapter.h"

namespace Asn::Math
{
    template <typename Space>
    class KohnShamAdapter;

    template <typename Space>
    class KohnShamAdapter final : public Adapter<Space>
    {
    protected:
        Matrix<Real> poisson_solution;
        Tensor<Real, 2> poisson_jump_term;
        Tensor<Real, 2> poisson_rhs;

        List<Tensor<Real, 2>> rho;
        Matrix<Real> kohn_sham_solution;
        Tensor<Real, 2> kohn_sham_jump_term;
        Vector<Real> energy;
        List<Real> num_elec;
        Tensor<Real, 2> ext_value;
        Tensor<Real, 2> har_value;
        Tensor<Real, 2> vxc_value;

    public:
        explicit KohnShamAdapter(Space &space);
        KohnShamAdapter(KohnShamAdapter &adapter) = delete;
        virtual ~KohnShamAdapter() override = default;

        void set_poisson_rhs(const Tensor<Real, 2> &rhs);
        void set_poisson_solution(const Matrix<Real> &solution);
        void set_kohn_sham_solution(const Matrix<Real> &solution, const List<Real> &num_elec);
        void set_energy(const Vector<Real> &energy);
        void set_ext_value(const Tensor<Real, 2> &ext_value);
        void set_har_value(const Tensor<Real, 2> &har_value);
        void set_vxc_value(const Tensor<Real, 2> &vxc_value);

        virtual Real h_indicator(const Int &pid) const override;

        KohnShamAdapter operator=(KohnShamAdapter &adapter) = delete;
    };
};

#endif
