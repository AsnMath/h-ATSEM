#include <iostream>

#include "../../include/Solver/Solver.hpp"
#include "../../include/Space/Space.hpp"
#include "../../include/Space/Instance/EllipticUtils.hpp"
#include "../../include/Space/Instance/EllipticAdapter.hpp"

using namespace Asn;
using namespace Asn::Math;

inline Real u(const Point<3> &x);
inline Real u_b(const Point<3> &x);
inline Real f(const Point<3> &x);

const String mesh_file = "../mesh/cube-1.mesh"; // Initial mesh.
const String out_mesh_file = "output.mesh"; // Saved mesh.
constexpr Int max_order = 4; // Degree of basis polynomial.
constexpr Int quad_order = 24; // Order of quadrature.
constexpr Real pcg_tol = 1e-15; // Tolerence of PCG.
constexpr Int max_iter = 10; // Max iteration of adaptive.
const String h_mode = "rate_max"; // h-adaptive mode.
constexpr Real h_rate = 0.5; // h-adaptive rate.
constexpr Int num_thread = -1; // Number of thread, less than or equal to 0 to use all cores.

int main()
{
    srand(0);

    if constexpr (num_thread > 0)
    {
        omp_set_num_threads(num_thread);
        Eigen::setNbThreads(num_thread);
        std::cout << "Num of thread is set to " << num_thread << std::endl;
    }

    std::cout << "Options:\n";
    std::cout << "    Initial Mesh: " << mesh_file << ".\n";
    std::cout << "    Degree of Basis Polynomial: " << max_order << std::endl;
    std::cout << "    Order of Quadrature: " << quad_order << std::endl;
    std::cout << "    Tolerence of PCG: " << pcg_tol << std::endl;
    std::cout << "    Max Iteration of Adaptive: " << max_iter << std::endl;
    std::cout << "    h-Adaptive Mode: " << h_mode << std::endl;
    std::cout << "    h-Adaptive Rate: " << h_rate << std::endl;
    std::cout << "    Number of Thread: " << num_thread << std::endl;
    std::cout << std::endl;

    Space<> space(mesh_file, max_order, quad_order);
    EllipticUtils<Space<>> utils(space);
    EllipticAdapter<Space<>> adapter(space);

    PCG pcg;
    PointJacobiPreconditioner preconditioner;

    for (Int iter = 1; iter <= max_iter; iter++)
    {
        std::cout << std::endl;
        std::cout << "Iteration " << iter << std::endl;
        std::cout << "    Building linear system ... ";
        utils.init();
        std::cout << "n_dof = " << space.num_dof << ", n_poly = " << space.num_poly() << std::endl;
        Tensor<Real, 2> f_value = space.get_value(f);
        SparseMatrix<Real> stiff_matrix = utils.get_stiff_matrix(MatrixMode::Full);
        Vector<Real> boundary_vec = utils.get_boundary(u_b);
        const Vector<Real> rhs = utils.get_rhs(f_value);
        utils.apply_boundary(stiff_matrix, boundary_vec);

        std::cout << "    Solving equation ... ";
        preconditioner.set_lhs(&stiff_matrix);
        pcg.set_rel_tol(pcg_tol);
        pcg.set_preconditioner(&preconditioner);
        const Vector<Real> solution = pcg.solve(stiff_matrix, boundary_vec + rhs, Vector<Real>::Zero(space.num_dof));
        switch (pcg.get_exit_type())
        {
            case RelTol:
                std::cout << " PCG reaches relative error with " << pcg.get_num_iter() << " iterations.\n";
                break;
            case AbsTol:
                std::cout << " PCG reaches Absolute error with " << pcg.get_num_iter() << " iterations.\n";
                break;
            default:
                std::cout << "            PCG reaches max iteration\n";
                break;
        }

        std::cout << "    Computing error ... ";
        const Tensor<Real, 2> tmp = space.get_value(solution) - space.get_value(u);
        const Real error = std::sqrt(space.inner_product(tmp, tmp));
        std::cout << "Error = " << error << std::endl;
        std::cout << "    Applying h-adaptive ... ";
        adapter.set_solution(solution);
        adapter.set_rhs(f_value);
        adapter.set_h_adaptive_flag(h_mode, h_rate);
        space.h_adaptive();
        space.init();
        space.check_regular();
    }
    space.save(out_mesh_file);
    return 0;
}

constexpr Real boundary_tol = 1e-7;

constexpr Real EXP_R_K = -10.0;

template <typename Type>
inline Type temp_u(const Point<3, Type> &x)
{
    const Type r = x.template lpNorm<2>();
    return exp(EXP_R_K * r);
}

inline Real u(const Point<3> &x)
{
    return temp_u<Real>(x);
}

inline Real u_b(const Point<3> &x)
{
    assert(x(0) <= -1.0 + boundary_tol || x(0) >= 1.0 - boundary_tol ||
            x(1) <= -1.0 + boundary_tol || x(1) >= 1.0 - boundary_tol ||
            x(2) <= -1.0 + boundary_tol || x(2) >= 1.0 - boundary_tol);
    return u(x);
}

inline Real f(const Point<3> &x)
{
    const Real r = x.lpNorm<2>();
    return -EXP_R_K * (2.0 + EXP_R_K * r) * exp(EXP_R_K * r) / r;
}
