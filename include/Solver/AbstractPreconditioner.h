#ifndef ASN_SOLVER_ABSTRACTPRECONDITIONER_H
#define ASN_SOLVER_ABSTRACTPRECONDITIONER_H

#include "../Config.h"

namespace Asn::Math
{
    class AbstractPreconditioner;

    class AbstractPreconditioner
    {
    public:
        AbstractPreconditioner() = default;
        AbstractPreconditioner(const AbstractPreconditioner &p) = default;
        virtual ~AbstractPreconditioner() = default;

        virtual AbstractPreconditioner &set_lhs(const SparseMatrix<Real> *lhs, const Real &tol = DEFAULT_ABS_TOL) = 0;
        virtual AbstractPreconditioner &set_rhs(const SparseMatrix<Real> *rhs, const Real &tol = DEFAULT_ABS_TOL) = 0;
        virtual AbstractPreconditioner &set_sigma(const Vector<Real> &sigma) = 0;

        virtual Matrix<Real> solve(const Matrix<Real> &x) const = 0;

        AbstractPreconditioner &operator=(const AbstractPreconditioner &p) = default;
    };
};

#endif


