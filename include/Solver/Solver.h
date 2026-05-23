#ifndef ASN_SOLVER_SOLVER_H
#define ASN_SOLVER_SOLVER_H

#include "AbstractSolver.h"

namespace Asn::Math
{
    class PCG;
    class LOBPCG;

    class PCG final : public AbstractSolver
    {
    public:
        PCG() = default;
        PCG(const PCG &solver) = default;
        virtual ~PCG() override = default;

        Vector<Real> solve(const SparseMatrix<Real> &A, const Vector<Real> &b);
        Vector<Real> solve(const SparseMatrix<Real> &A, const Vector<Real> &b, const Vector<Real> &init_guess);

        PCG &operator=(const PCG &solver) = default;
    };

    class LOBPCG final : public AbstractSolver
    {
    public:
        LOBPCG() = default;
        LOBPCG(const LOBPCG &solver) = default;
        virtual ~LOBPCG() override = default;

        Pair<Vector<Real>, Matrix<Real>> solve(const SparseMatrix<Real> &A, const SparseMatrix<Real> &B, const Matrix<Real> &init_guess);

        LOBPCG &operator=(const LOBPCG &solver) = default;
    };
};

#endif
