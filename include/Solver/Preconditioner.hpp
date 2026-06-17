#ifndef ASN_SOLVER_PRECONDITIONER_HPP
#define ASN_SOLVER_PRECONDITIONER_HPP

#include "Preconditioner.h"

namespace Asn::Math
{
    inline PointJacobiPreconditioner::PointJacobiPreconditioner()
    {
        this->init_rhs = false;
        this->init_lhs = false;
        this->init_sigma = false;
    }

    inline PointJacobiPreconditioner &PointJacobiPreconditioner::set_lhs(const SparseMatrix<Real> *lhs, const Real &tol)
    {
        assert(tol >= 0.0);
        this->lhs = lhs->diagonal();
        for (Int i = 0; i < this->lhs.rows(); ++i)
        {
            if (std::abs(this->lhs(i)) < tol)
            {
                this->lhs(i) = 1.0;
            }
        }
        this->init_lhs = true;
        return *this;
    }

    inline PointJacobiPreconditioner &PointJacobiPreconditioner::set_rhs(const SparseMatrix<Real> *rhs, const Real &tol)
    {
        assert(tol >= 0.0);
        this->rhs = rhs->diagonal();
        for (Int i = 0; i < this->rhs.rows(); ++i)
        {
            if (std::abs(this->rhs(i)) < tol)
            {
                this->rhs(i) = 1.0;
            }
        }
        this->init_rhs = true;
        return *this;
    }

    inline PointJacobiPreconditioner &PointJacobiPreconditioner::set_sigma(const Vector<Real> &sigma)
    {
        this->sigma = sigma;
        this->init_sigma = true;
        return *this;
    }

    inline Matrix<Real> PointJacobiPreconditioner::solve(const Matrix<Real> &x) const
    {
        if (this->init_lhs)
        {
            if (this->init_rhs && this->init_sigma)
            {
                Matrix<Real> res(x.rows(), x.cols());
                for (Int i = 0; i < x.cols(); i++)
                {
                    res.col(i) = (this->lhs - this->sigma(i) * this->rhs).cwiseInverse().asDiagonal() * x.col(i);
                }
                return res;
            }
            else
            {
                return this->lhs.cwiseInverse().asDiagonal() * x;
            }
        }
        ASN_ERROR("The preconditioner is not initialized.");
    }
}

#endif

