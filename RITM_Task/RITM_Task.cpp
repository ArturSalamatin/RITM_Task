// RITM_Task.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <tuple>

#include "ConcreteSolver.h"
#include "LocalSolver/RK2Explicit.h"
#include "LocalSolver/RK4Explicit.h"
#include "OutputStream/OutputStream.h"

template<typename rk>
 auto //std::pair<ConcreteVerifier, ConcreteSolver>
    verify(
        const ConcreteVerifier& verifier_temp,
        const ConcreteSolver& solver_temp, 
        const ConcreteParams& params)
{
    auto solver{ solver_temp };
    solver.solve(rk{params});
    auto verifier{ verifier_temp };
    verifier.verify(solver);
    return std::pair{ verifier, solver };
}

int main()
{
    double C{ 1.0 };
    double L{ 1.0 };
    double omega2{ C / L };
    double U0{ 1.0 };
    double I0{ 0.0 };
    double q0{ U0 / C };

    ConcreteParams params{ omega2 };
    ConcreteState init_state{U0/C, I0};
    ConcreteRHS rhs{};
    // for test purposes
    // auto init_differential{ rhs(init_state, params) };
    UniformGrid grid{0.0, 100.0, 0.1};
    ConcreteSolver solver{ params, grid, init_state, rhs };
    ConcreteVerifier verifier{ C, L, q0, I0 };
    auto [verifier4, solver4] = 
        verify<CauchySolver::LocalSolver::RK4<ConcreteState, ConcreteParams>>(verifier, solver, params);
    auto [verifier2, solver2] {
        verify<CauchySolver::LocalSolver::RK2<ConcreteState, ConcreteParams>>(verifier, solver, params)};

    CauchySolver::Printer 
        printer4{"output4.txt"},
        printer2{"output2.txt"};

    printer4.print(verifier4, solver4);
    printer2.print(verifier2, solver2);
}
