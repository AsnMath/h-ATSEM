#include <iostream>

#include "../../include/Solver/Solver.hpp"
#include "../../include/Space/Space.hpp"
#include "../../include/Space/Spline.hpp"
#include "../../include/Space/Instance/KohnShamUtils.hpp"
#include "../../include/Space/Instance/KohnShamAdapter.hpp"
#include "../../include/Space/Instance/EllipticAdapter.hpp"

using namespace Asn;
using namespace Asn::Math;

const String mesh_file = "../mesh/cube-10.mesh"; // Initial mesh.
const String out_mesh_file = "output.mesh"; // Saved mesh.
constexpr Int max_order = 4; // Degree of basis polynomial.
constexpr Int quad_order = 24; // Order of quadrature.
constexpr Real pcg_tol = 1e-8; // Tolerence of PCG.
constexpr Real lobpcg_tol = 1e-6; // Tolerence of LOBPCG.
constexpr Real scf_tol = 1e-4; // Tolerence of SCF iteration.
constexpr Int max_iter = 8; // Max iteration of adaptive.
const String h_mode = "rate_max"; // h-adaptive mode.
constexpr Real h_rate = 0.5; // h-adaptive rate.
constexpr Int num_thread = -1; // Number of thread, less than or equal to 0 to use all cores.
constexpr Real factor_mixture = 0.5; // Parameter for Pulay's mixture.
constexpr Int num_history = 5; // Parameter for Pulay's mixture.
constexpr Bool use_interp = true; // Using interpolation.

const List<Point<3>> nuclei = {{0.0, 0.0, 0.0}};
const List<Real> nuclei_charge = {2.0};
const List<Real> num_elec = {2.0};
const List<Real> num_elec_up = {1.0};
const List<Real> num_elec_down = {1.0};
constexpr Real real_energy = -2.834836;

int main()
{
    srand(0);

    if (num_thread > 0)
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
    std::cout << "    Tolerence of LOBPCG: " << lobpcg_tol << std::endl;
    std::cout << "    Tolerence of SCF Iteration: " << scf_tol << std::endl;
    std::cout << "    Max Iteration of Adaptive: " << max_iter << std::endl;
    std::cout << "    h-Adaptive Mode: " << h_mode << std::endl;
    std::cout << "    h-Adaptive Rate: " << h_rate << std::endl;
    std::cout << "    Number of Thread: " << num_thread << std::endl;
    std::cout << "    Using Interpolation: " << use_interp << std::endl;
    std::cout << std::endl;

    Space<> space(mesh_file, max_order, quad_order);
    KohnShamUtils<Space<>> utils(space);
    KohnShamAdapter<Space<>> adapter(space);
    EllipticAdapter<Space<>> poisson_adapter(space);

    Math::PointJacobiPreconditioner pcg_preconditioner;
    Math::PointJacobiPreconditioner lobpcg_preconditioner;
    Math::PCG pcg;
    Math::LOBPCG lobpcg;
    pcg.set_rel_tol(pcg_tol).set_preconditioner(&pcg_preconditioner);
    lobpcg.set_rel_tol(lobpcg_tol).set_max_iter(1024).set_preconditioner(&lobpcg_preconditioner);

    Spline spline;
    Matrix<Real> prev_solution;

    const Real E_ION_ION = utils.get_energy_ion_ion(nuclei, nuclei_charge);

    for (Int iter = 1; iter <= max_iter; iter++)
    {
        utils.init();

        std::cout << std::endl;
        std::cout << "Iteration " << iter << " n_dof = " << space.num_dof << ", n_poly = " << space.num_poly() << std::endl;

        Vector<Real> energy = Vector<Real>::Ones(num_elec.size());
        Vector<Real> epsilon = Vector<Real>::Ones(num_elec.size());
        Matrix<Real> solution;
        Vector<Real> har_solution = Vector<Real>::Random(space.num_dof);

        Tensor<Real, 2> rho;

        if (iter == 1 || (!use_interp))
        {
            solution = Matrix<Real>::Random(space.num_dof, num_elec.size());
            rho = utils.get_rho([](const Real &r) { return std::exp(-2.0 * r) / std::numbers::pi; }, nuclei, nuclei_charge);
        }
        else
        {
            solution = spline.interp(prev_solution, space);
            rho = utils.get_rho(solution, num_elec);
        }

        Tensor<Real, 2> ext = utils.get_ext(nuclei, nuclei_charge);
        Tensor<Real, 2> har;
        Tensor<Real, 2> exc;
        Tensor<Real, 2> vxc;

        SparseMatrix<Real> ks_stiff_matrix = utils.get_stiff_matrix(MatrixMode::DelBoundary) / 2.0;
        SparseMatrix<Real> ks_mass_matrix = utils.get_mass_matrix(MatrixMode::DelBoundary);

        SparseMatrix<Real> stiff_matrix = utils.get_stiff_matrix(MatrixMode::Full);
        Vector<Real> boundary_vec = Vector<Real>::Zero(stiff_matrix.rows());
        utils.apply_boundary(stiff_matrix, boundary_vec);

        pcg_preconditioner.set_lhs(&stiff_matrix);

        lobpcg_preconditioner.set_lhs(&ks_stiff_matrix);
        lobpcg_preconditioner.set_rhs(&ks_mass_matrix);

        List<Tensor<Real, 2>> rho_history;
        List<Tensor<Real, 2>> res_history;

        Int scf_iter = 0;
        while (true)
        {
            scf_iter++;
            std::cout << "    SCF Iteration " << scf_iter << std::endl;
            std::cout << "        Solving Hartree ...";
            stiff_matrix = utils.get_stiff_matrix(MatrixMode::Full);
            auto [q, p, Q, O] = utils.get_multipole_expansion(rho);
            boundary_vec = utils.get_boundary([&p, &q, &Q, &O](const Vector<Real> &x)
            {
                const Real r = x.lpNorm<2>();
                const Real r3 = r * r * r;
                const Real r5 = r3 * r * r;
                const Real r7 = r5 * r * r;
                Real res = q / r;
                res += p.dot(x) / (r3);
                res += x.dot(Q * x) / (2.0 * r5);
                for (Int i = 0; i < 3; i++)
                {
                    for (Int j = 0; j < 3; j++)
                    {
                        for (Int k = 0; k < 3; k++)
                        {
                            res += ((O(i, j, k) * x(i) * x(j) * x(k)) / (6.0 * r7));
                        }
                    }
                }
                return res;
            });
            utils.apply_boundary(stiff_matrix, boundary_vec);
            const Vector<Real> rhs = utils.get_rhs(4.0 * std::numbers::pi * rho);
            har_solution = pcg.solve(stiff_matrix, boundary_vec + rhs, har_solution);
            switch (pcg.get_exit_type())
            {
                case RelTol:
                    std::cout << "PCG reaches relative error with " << pcg.get_num_iter() << " iterations.\n";
                    break;
                case AbsTol:
                    std::cout << "PCG reaches Absolute error with " << pcg.get_num_iter() << " iterations.\n";
                    break;
                default:
                    std::cout << "PCG reaches max iteration\n";
                    break;
            }

            har = space.get_value(har_solution);
            std::tie(vxc, exc) = utils.get_lda(rho);
            const SparseMatrix<Real> A = ks_stiff_matrix + utils.get_effect_matrix(ext + har + vxc);

            const Vector<Real> new_energy = utils.get_energy(epsilon, num_elec, solution, har, vxc, exc);
            const Real scf_energy_error = (new_energy - energy).lpNorm<2>();
            energy = new_energy;
            if (scf_energy_error <= scf_tol)
            {
                break;
            }

            std::cout << "        Current energy is " << energy << std::endl;
            std::cout << "        Solving Kohn-Sham ...";
            auto [eps, sol] = lobpcg.solve(A, ks_mass_matrix, space.map_matrix_row(solution, utils.dof_index));
            switch (lobpcg.get_exit_type())
            {
                case RelTol:
                    std::cout << " LOBPCG reaches relative error with " << lobpcg.get_num_iter() << " iterations.\n";
                    break;
                case AbsTol:
                    std::cout << " LOBPCG reaches Absolute error with " << lobpcg.get_num_iter() << " iterations.\n";
                    break;
                default:
                    std::cout << " LOBPCG reaches max iteration\n";
                    break;
            }

            epsilon = eps;
            solution = space.map_matrix_row(sol, utils.index_dof, space.num_dof);

            const Tensor<Real, 2> new_rho = utils.get_rho(solution, num_elec);
            rho_history.push_back(rho);
            res_history.push_back(new_rho - rho);
            if (static_cast<Int>(rho_history.size()) > num_history)
            {
                rho_history.erase(rho_history.begin());
                res_history.erase(res_history.begin());
            }
            rho = utils.pulay_mixing(rho_history, res_history, factor_mixture);
        }

        const Real total_energy = energy.sum() + E_ION_ION;

        std::cout << "    Total energy = " << to_fixed_precision(total_energy, 6) << " / " << to_fixed_precision(real_energy, 6) <<
                " with " << scf_iter << " SCF Iterations" << std::endl;

        if (use_interp)
        {
            spline.init(space);
            prev_solution = solution;
        }
        std::cout << "    Applying h-adaptive ... \n";
        space.clear_h_adaptive_flag();

        poisson_adapter.set_solution(har_solution);
        poisson_adapter.set_rhs(4.0 * std::numbers::pi * rho);
        poisson_adapter.set_h_adaptive_flag(h_mode, h_rate);

        adapter.set_poisson_rhs(4.0 * std::numbers::pi * rho);
        adapter.set_poisson_solution(har_solution);
        adapter.set_kohn_sham_solution(solution, num_elec);
        adapter.set_energy(energy);
        adapter.set_ext_value(ext);
        adapter.set_har_value(har);
        adapter.set_vxc_value(vxc);

        adapter.set_h_adaptive_flag(h_mode, h_rate);
        space.h_adaptive();
        space.init();

        space.check_regular();
    }
    space.save(out_mesh_file);
    return 0;
}
