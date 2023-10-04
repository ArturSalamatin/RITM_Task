#pragma once
#include <vector>

#include "AbstractSolver/AbstractSolver.h"
#include "AbstractSolver/State.h"
#include "Grids/UniformGrid.h"

// The state is described by two variables
// y = (q; I) --- (Charge; Current)
constexpr size_t problem_size = 2;
using ConcreteState = 
	CauchySolver::AbstractSolver::State<problem_size>;

// Parameters, descriing 
// the rhs of the problem
struct ConcreteParams
{
	double omega2;
};

// The rhs f(y) describes the 
// dynamics of the system,
// the evolution 
// of charge, f[0], and 
// of current, f[1]
struct ConcreteRHS
{
	ConcreteState operator()(
		const ConcreteState& y,
		const ConcreteParams& params) const
	{
		return ConcreteState
		{
			y[1],
			-params.omega2 * y[0]
		};
	}
};

// Define a solver for the particular problem.
// The local solver, i.e., RK4 ir RK2,
// are not incorporated into the solver yet.
using ConcreteSolver =
CauchySolver::AbstractSolver::AbstractSolver
<
	ConcreteState,
	ConcreteParams,
	UniformGrid,
	ConcreteRHS
>;




// Verify the numerical solution against the 
// analytical solution.
struct ConcreteVerifier
{
	ConcreteVerifier(
		double C, double L,
		double q0, double I0) :
		C{ C }, L{ L }
	{
		E = (C * q0 * q0 + L * I0 * I0) / 2.0;
	}

	void verify(const ConcreteSolver& solver)
	{
		calc_energy.reserve(solver.grid.size());
		abs_diff.reserve(solver.grid.size());
		rel_diff.reserve(solver.grid.size());
		for (size_t i = 0; i < calc_energy.capacity(); ++i)
		{
			double q{ solver.result(i)[0] };
			double I{ solver.result(i)[1] };
			calc_energy.push_back((C * q * q + L * I * I) / 2.0);
			abs_diff.push_back(std::abs(calc_energy.back() - E));
			rel_diff.push_back(abs_diff.back() / std::abs(calc_energy.back() + E));
		}
	}

	double C, L;
	double E;

	std::vector<double> calc_energy;
	std::vector<double> abs_diff, rel_diff;
};

