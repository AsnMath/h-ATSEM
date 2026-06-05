#ifndef ASN_SOLVER_SOLVER_HPP
#define ASN_SOLVER_SOLVER_HPP

#include "Solver.h"

#include "AbstractSolver.hpp"
#include "Preconditioner.hpp"

namespace Asn::Math
{
    constexpr Real MAX_MODIFY_TOL = 1e16;

    inline Matrix<Real> safe_decomposition(const Matrix<Real> &A)
    {
        const Int N = A.rows();
        Real sigma = std::numeric_limits<Real>::epsilon() * A.norm();
        Matrix<Real> O = A;
        while (sigma <= MAX_MODIFY_TOL)
        {
            if (const Eigen::LDLT<Matrix<Real>> chol(O); chol.info() == Eigen::ComputationInfo::Success)
            {
                const Matrix<Real> U = chol.matrixU().solve(Matrix<Real>::Identity(N, N));
                if (const Vector<Real> D = chol.vectorD(); D.minCoeff() > 0.0)
                {
                    return chol.transpositionsP().transpose() * U * D.cwiseSqrt().cwiseInverse().asDiagonal();
                }
            }
            sigma *= 10.0;
            O.diagonal().array() += sigma;
        }
        ASN_ERROR("Function safe_decomposition fails");
    }

    inline Tuple<Vector<Real>, Matrix<Real>> safe_eigen_solver(const Matrix<Real> &A, const Matrix<Real> &B)
    {
        const Int N = A.rows();
        Real sigma = std::numeric_limits<Real>::epsilon() * B.norm();
        Matrix<Real> O = B;
        while (sigma <= MAX_MODIFY_TOL)
        {
            if (const Eigen::GeneralizedSelfAdjointEigenSolver<Matrix<Real>> eigs(A, O); eigs.info() == Eigen::ComputationInfo::Success)
            {
                return std::make_tuple(eigs.eigenvalues(), eigs.eigenvectors());
            }
            sigma *= 10.0;
            O.diagonal().array() += sigma;
        }
        ASN_ERROR("Function safe_eigen_solver fails");
    }

    inline Matrix<Real> stack_4_matricies(const Matrix<Real> &A11, const Matrix<Real> &A12, const Matrix<Real> &A22)
    {
        Matrix<Real> res(A11.rows() + A22.rows(), A11.cols() + A22.cols());
        res.topLeftCorner(A11.rows(), A11.cols()) = A11;
        res.topRightCorner(A12.rows(), A12.cols()) = A12;
        res.bottomRightCorner(A22.rows(), A22.cols()) = A22;
        res = res.selfadjointView<Eigen::Upper>();
        return res;
    }

    inline Matrix<Real> stack_9_matricies(const Matrix<Real> &A11, const Matrix<Real> &A12, const Matrix<Real> &A13,
                                          const Matrix<Real> &A22, const Matrix<Real> &A23, const Matrix<Real> &A33)
    {
        Matrix<Real> res(A11.rows() + A22.rows() + A33.rows(), A11.cols() + A22.cols() + A33.cols());
        res.block(0, 0, A11.rows(), A11.cols()) = A11;
        res.block(0, A11.cols(), A12.rows(), A12.cols()) = A12;
        res.block(0, A11.cols() + A12.cols(), A13.rows(), A13.cols()) = A13;
        res.block(A11.rows(), A11.cols(), A22.rows(), A22.cols()) = A22;
        res.block(A11.rows(), A11.cols() + A12.cols(), A23.rows(), A23.cols()) = A23;
        res.block(A11.rows() + A22.rows(), A11.cols() + A12.cols(), A33.rows(), A33.cols()) = A33;
        res = res.selfadjointView<Eigen::Upper>();
        return res;
    }

    inline void sort_epairs(Vector<Real> &evalues, Matrix<Real> &evectors)
    {
        MultiMap<Real, Vector<Real>> epairs;
        for (Int i = 0; i < evectors.cols(); i++)
        {
            epairs.insert(std::make_pair(evalues(i), evectors.col(i)));
        }
        Int i = 0;
        for (const auto &[fst, snd] : epairs)
        {
            evalues(i) = fst;
            evectors.col(i) = snd;
            i++;
        }
        return;
    }

    inline Vector<Real> PCG::solve(const SparseMatrix<Real> &A, const Vector<Real> &b)
    {
        return solve(A, b, Vector<Real>::Zero(b.rows()));
    }

    inline Vector<Real> PCG::solve(const SparseMatrix<Real> &A, const Vector<Real> &b, const Vector<Real> &init_guess)
    {
        const Int N = b.rows();
        const Real r_tol = rel_tol * b.lpNorm<2>();
        Vector<Real> w(N), x(N), r(N), z(N), p(N);

        x = init_guess;
        r = b - A * x;
        if (precond != nullptr)
        {
            z.noalias() = precond->solve(r);
        }
        else
        {
            z.noalias() = r;
        }
        p = z;
        Real rho = r.dot(z);
        for (num_iter = 0; num_iter < max_iter; num_iter++)
        {
            const Real norm = r.norm();
            if (norm <= r_tol)
            {
                exit_type = RelTol;
                break;
            }
            if (norm <= abs_tol)
            {
                exit_type = AbsTol;
                break;
            }
            w.noalias() = A * p;
            const Real alpha = rho / p.dot(w);
            x.noalias() += alpha * p;
            r.noalias() -= alpha * w;
            if (precond != nullptr)
            {
                z.noalias() = precond->solve(r);
            }
            else
            {
                z.noalias() = r;
            }
            const Real t_rho = rho;
            rho = r.dot(z);
            const Real beta = rho / t_rho;
            p = z + beta * p;
        }
        if (num_iter >= max_iter)
        {
            exit_type = MaxIter;
        }
        return x;
    }

    inline Pair<Vector<Real>, Matrix<Real>> LOBPCG::solve(const SparseMatrix<Real> &A, const SparseMatrix<Real> &B, const Matrix<Real> &init_guess)
    {
        const Int M = init_guess.cols();
        const Real norm_a = A.norm();
        const Real norm_b = B.norm();

        Matrix<Real> X = init_guess;
        Matrix<Real> W, AW, BW, P, AP, BP, AX, BX, inv_R, C;
        Vector<Real> lambda(M);

        BX = B * X;
        inv_R = safe_decomposition(X.transpose() * BX);
        X = X * inv_R;
        BX = BX * inv_R;
        AX = A * X;
        Eigen::SelfAdjointEigenSolver<Matrix<Real>> eigs(Matrix<Real>(X.transpose() * AX));
        assert(eigs.info() == Eigen::ComputationInfo::Success);
        lambda = eigs.eigenvalues();
        Matrix<Real> TMP = eigs.eigenvectors();
        sort_epairs(lambda, TMP);
        X = X * TMP;
        AX = AX * TMP;
        BX = BX * TMP;

        for (num_iter = 0; num_iter < max_iter; num_iter++)
        {
            W = AX - BX * lambda.asDiagonal();
            Bool end_rel = true;
            Bool end_abs = true;
            for (Int i = 0; i < M; i++)
            {
                const Real norm = W.col(i).lpNorm<2>();
                const Bool fail_rel = norm > rel_tol * (norm_a + std::abs(lambda(i)) * norm_b);
                const Bool fail_abs = norm > abs_tol;
                if (fail_rel)
                {
                    end_rel = false;
                }
                if (fail_abs)
                {
                    end_abs = false;
                }
                if (fail_rel && fail_abs)
                {
                    break;
                }
            }
            if (end_rel)
            {
                exit_type = RelTol;
                break;
            }
            if (end_abs)
            {
                exit_type = AbsTol;
                break;
            }
            if (precond != nullptr)
            {
                precond->set_sigma(lambda);
                W = precond->solve(W);
            }
            BW = B * W;
            inv_R = safe_decomposition(W.transpose() * BW);
            W = W * inv_R;
            BW = BW * inv_R;
            AW = A * W;
            Matrix<Real> gram_A, gram_B;
            if (num_iter > 0)
            {
                inv_R = safe_decomposition(P.transpose() * BP);
                P = P * inv_R;
                AP = AP * inv_R;
                BP = BP * inv_R;
                gram_A = stack_9_matricies(X.transpose() * AX, X.transpose() * AW, X.transpose() * AP,
                                           W.transpose() * AW, W.transpose() * AP,
                                           P.transpose() * AP);
                gram_B = stack_9_matricies(X.transpose() * BX, X.transpose() * BW, X.transpose() * BP,
                                           W.transpose() * BW, W.transpose() * BP,
                                           P.transpose() * BP);
            }
            else
            {
                gram_A = stack_4_matricies(X.transpose() * AX, X.transpose() * AW,
                                           W.transpose() * AW);
                gram_B = stack_4_matricies(X.transpose() * BX, X.transpose() * BW,
                                           W.transpose() * BW);
            }
            std::tie(lambda, C) = safe_eigen_solver(gram_A, gram_B);
            sort_epairs(lambda, C);
            lambda = lambda.topRows(M).eval();
            C = C.leftCols(M).eval();
            const Matrix<Real> CX = C.block(0, 0, M, M);
            const Matrix<Real> CW = C.block(M, 0, M, M);
            if (num_iter > 0)
            {
                const Matrix<Real> CP = C.block(M + M, 0, M, M);
                P = W * CW + P * CP;
                AP = AW * CW + AP * CP;
                BP = BW * CW + BP * CP;
            }
            else
            {
                P = W * CW;
                AP = AW * CW;
                BP = BW * CW;
            }
            X = X * CX + P;
            AX = AX * CX + AP;
            BX = BX * CX + BP;
        }
        if (num_iter >= max_iter)
        {
            exit_type = MaxIter;
        }
        sort_epairs(lambda, X);
        return std::make_pair(lambda, X);
    }
};

#endif
