#ifndef ASN_SOLVER_ABSTRACTSOLVER_HPP
#define ASN_SOLVER_ABSTRACTSOLVER_HPP

#include "AbstractSolver.h"

namespace Asn::Math
{
    inline AbstractSolver::AbstractSolver() : num_iter(0), precond(nullptr)
    {
        this->abs_tol = DEFAULT_ABS_TOL;
        this->rel_tol = DEFAULT_REL_TOL;
        this->max_iter = DEFAULT_MAX_ITER;
        this->exit_type = Unknown;
    }

    inline AbstractSolver &AbstractSolver::set_abs_tol(const Real &abs_tol)
    {
        assert(abs_tol >= 0.0);
        this->abs_tol = abs_tol;
        return *this;
    }

    inline AbstractSolver &AbstractSolver::set_rel_tol(const Real &rel_tol)
    {
        assert(rel_tol >= 0.0);
        this->rel_tol = rel_tol;
        return *this;
    }

    inline AbstractSolver &AbstractSolver::set_max_iter(const Int &max_iter)
    {
        assert(max_iter >= 0.0);
        this->max_iter = max_iter;
        return *this;
    }

    inline AbstractSolver &AbstractSolver::set_preconditioner(AbstractPreconditioner *precond)
    {
        this->precond = precond;
        return *this;
    }

    inline ExitType AbstractSolver::get_exit_type() const
    {
        return this->exit_type;
    }

    inline Int AbstractSolver::get_num_iter() const
    {
        return this->num_iter;
    }
};

#endif
