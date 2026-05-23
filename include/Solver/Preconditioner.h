#ifndef ASN_SOLVER_PRECONDITIONER_H
#define ASN_SOLVER_PRECONDITIONER_H

#include "AbstractPreconditioner.h"

namespace Asn::Math
{
    class PointJacobiPreconditioner;

    class PointJacobiPreconditioner final : public AbstractPreconditioner
    {
    private:
        Bool init_lhs;
        Bool init_rhs;
        Bool init_sigma;

        Vector<Real> lhs;
        Vector<Real> rhs;
        Vector<Real> sigma;

    public:
        PointJacobiPreconditioner();
        PointJacobiPreconditioner(const PointJacobiPreconditioner &p) = default;
        virtual ~PointJacobiPreconditioner() override = default;

        virtual PointJacobiPreconditioner &set_lhs(const SparseMatrix<Real> *lhs, const Real &tol = DEFAULT_ABS_TOL) override;
        virtual PointJacobiPreconditioner &set_rhs(const SparseMatrix<Real> *rhs, const Real &tol = DEFAULT_ABS_TOL) override;
        virtual PointJacobiPreconditioner &set_sigma(const Vector<Real> &sigma) override;

        virtual Matrix<Real> solve(const Matrix<Real> &x) const override;

        PointJacobiPreconditioner &operator=(const PointJacobiPreconditioner &p) = default;
    };
};

#endif
