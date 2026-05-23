#ifndef ASN_SOLVER_ABSTRACTSOLVER_H
#define ASN_SOLVER_ABSTRACTSOLVER_H

#include "AbstractPreconditioner.h"

namespace Asn::Math
{
    enum ExitType
    {
        Success,
        Unknown,
        AbsTol,
        RelTol,
        MaxIter,
    };

    class AbstractSolver;

    class AbstractSolver
    {
    protected:
        Real abs_tol;
        Real rel_tol;
        Int max_iter;

        Int num_iter;
        ExitType exit_type;

        AbstractPreconditioner *precond;

    public:
        AbstractSolver();
        AbstractSolver(const AbstractSolver &solver) = default;
        virtual ~AbstractSolver() = default;

        AbstractSolver &set_abs_tol(const Real &abs_tol);
        AbstractSolver &set_rel_tol(const Real &rel_tol);
        AbstractSolver &set_max_iter(const Int &max_iter);
        AbstractSolver &set_preconditioner(AbstractPreconditioner *precond);

        ExitType get_exit_type() const;
        Int get_num_iter() const;

        AbstractSolver &operator=(const AbstractSolver &solver) = default;
    };
};

#endif
