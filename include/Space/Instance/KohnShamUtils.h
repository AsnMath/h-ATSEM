#ifndef ASN_MATH_SPACE_KOHN_SHAM_UTILS_H
#define ASN_MATH_SPACE_KOHN_SHAM_UTILS_H

#include "EllipticUtils.h"

namespace Asn::Math
{
    template <typename Space>
    class KohnShamUtils;

    template <typename Space>
    class KohnShamUtils final : public EllipticUtils<Space>
    {
    private:
        static constexpr Int LIBXC_LDA_X_ID = 1;
        static constexpr Int LIBXC_LDA_C_ID = 7;

    private:
        Real element_effect_matrix(const Int &pid, const Index<3> &i, const Index<3> &j, const Tensor<Real, 2> &eff) const;

    public:
        explicit KohnShamUtils(Space &space);
        KohnShamUtils(const KohnShamUtils &abstract_utils) = delete;
        virtual ~KohnShamUtils() override = default;

        virtual void init() override;

        Tensor<Real, 2> get_rho(const Func<Real(const Real &)> &rho, const List<Point<Space::MESH_DIM>> &pos, const List<Real> &charge) const;
        Tensor<Real, 2> get_rho(const Matrix<Real> &x, const List<Real> &F) const;

        Tensor<Real, 2> get_ext(const List<Point<Space::MESH_DIM>> &pos, const List<Real> &F) const;

        Tuple<Tensor<Real, 2>, Tensor<Real, 2>> get_lda(const Tensor<Real, 2> &rho) const;
        Tensor<Real, 2> get_lda_vxc(const Tensor<Real, 2> &rho) const;
        Tensor<Real, 2> get_lda_exc(const Tensor<Real, 2> &rho) const;

        Tuple<Tensor<Real, 2>, Tensor<Real, 2>, Tensor<Real, 2>> get_lsd(const Tensor<Real, 2> &rho_1, const Tensor<Real, 2> &rho_2) const;
        Pair<Tensor<Real, 2>, Tensor<Real, 2>> get_lsd_vxc(const Tensor<Real, 2> &rho_1, const Tensor<Real, 2> &rho_2) const;
        Tensor<Real, 2> get_lsd_exc(const Tensor<Real, 2> &rho_1, const Tensor<Real, 2> &rho_2) const;

        SparseMatrix<Real> get_effect_matrix(const Tensor<Real, 2> &eff) const;

        Tuple<Real, Vector<Real>, Matrix<Real>, Tensor<Real, 3>> get_multipole_expansion(const Tensor<Real, 2> &rho) const;

        static Real get_energy_ion_ion(const List<Point<Space::MESH_DIM>> &pos, const List<Real> &charge);

        Vector<Real> get_energy(const Vector<Real> &eps, const List<Real> &F, const Matrix<Real> &solution,
                                const Tensor<Real, 2> &har, const Tensor<Real, 2> &vxc, const Tensor<Real, 2> &exc) const;

        Tuple<Vector<Real>, Vector<Real>> get_energy(const Vector<Real> &eps_up, const Vector<Real> &eps_down,
                                                     const List<Real> &F_up, const List<Real> &F_down,
                                                     const Matrix<Real> &solution_up, const Matrix<Real> &solution_down, const Tensor<Real, 2> &har,
                                                     const Tensor<Real, 2> &vxc_up, const Tensor<Real, 2> &vxc_down, const Tensor<Real, 2> &exc) const;

        Tensor<Real, 2> pulay_mixing(const List<Tensor<Real, 2>> &rho, const List<Tensor<Real, 2>> &res, const Real &factor) const;

        KohnShamUtils &operator=(const KohnShamUtils &abstract_utils) = delete;
    };
};

#endif
