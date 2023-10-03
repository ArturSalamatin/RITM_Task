#pragma once
#include <array>
#include <vector>
#include <cmath>

namespace CauchySolver
{
	namespace AbstractSolver
	{
		template<size_t problem_size_t>
		struct State
		{
			template<typename... Args>
			State(Args... state) :
				state{ state... }
			{}

			State(const State& state) = default;

			State() = delete;

			auto operator[](size_t i) const
			{
				return state[i];
			}

			auto& operator[](size_t i)
			{
				return state[i];
			}

			std::array<double, problem_size_t>
				state;

			constexpr static size_t size()
			{
				return problem_size_t;
			}
		};

		template<size_t problem_size_t>
		State<problem_size_t> operator+(
			const State<problem_size_t>& lhs,
			const State<problem_size_t>& rhs)
		{
			State<problem_size_t> out{ lhs };
			for (size_t i = 0; i < out.size(); ++i)
				out[i] += rhs[i];
			return out;
		}

		template<size_t problem_size_t>
		State<problem_size_t> operator*(
			double lhs,
			const State<problem_size_t>& rhs)
		{
			State<problem_size_t> out{ rhs };
			for (size_t i = 0; i < out.size(); ++i)
				out[i] *= lhs;
			return out;
		}
		template<size_t problem_size_t>
		State<problem_size_t> operator*(
			const State<problem_size_t>& lhs,
			double rhs)
		{
			return rhs * lhs;
		}

		template<
			typename State_t,
			typename Params_t,
			typename Grid_t,
			typename RHS_t>
			struct Solver
		{
			Params_t params;
			Grid_t grid;
			State_t init_state;
			RHS_t rhs;

			Solver(
				const Params_t& params,
				const Grid_t& grid,
				const State_t& init_state,
				const RHS_t& rhs) :
				params{ params },
				grid{ grid },
				init_state{ init_state },
				rhs{ rhs }
			{
				its_result.reserve(grid.size());
			}

			template<
				//template<typename State_t, typename Params_t>
				typename LocalSolver_t>
				void solve(const LocalSolver_t& local_solver)
			{
				its_result.emplace_back(init_state);
				for (size_t i = 1; i < grid.size(); ++i)
					its_result.emplace_back(local_solver(its_result.back(), rhs, grid.step(i)));
			}

			const std::vector<State_t>& result() const
			{
				if (its_result.size() == 0)
					throw std::exception("One has to call solve() method to compute the result.");
				return its_result;
			}

			const State_t result(size_t i) const
			{
				return result()[i];
			}

		protected:
			std::vector<State_t> its_result;
		};
	} // AbstractSolver
} // CauchySolver

	struct UniformGrid
	{
		explicit UniformGrid(double start, double end, size_t node_count) :
			grid(node_count),
			its_step{(end-start)/(node_count-1ull)}
		{
			for (int i = 0; i < size(); ++i)
				grid[i] = start + step(i) * i;
		}

		explicit UniformGrid(double start, double end, double step) :
			UniformGrid{ 
			start, end, 
			static_cast<size_t>(std::ceil((end - start) / step) + 1ull) }
		{}

		UniformGrid(const UniformGrid&) = default;

		std::vector<double> grid;

		double start() const
		{
			return grid.front();
		}
		double end() const
		{
			return grid.back();
		}
		// general Grid interface must support the
		// case of non-uniform grid
		double step(size_t i) const
		{
			return its_step;
		}
		double size() const
		{
			return grid.size();
		}

		double its_step;
		
	};

	struct ConcreteParams
	{
		double omega2;
	};

	constexpr size_t problem_size = 2;
	using ConcreteState = CauchySolver::AbstractSolver::State<problem_size>;

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

	using ConcreteSolver =
		CauchySolver::AbstractSolver::Solver
		<
			ConcreteState,
			ConcreteParams,
			UniformGrid,
			ConcreteRHS
		>;

	template
		<
		typename State_t,
		typename Params_t
		>
		struct RK4
	{
		Params_t params;
		RK4(const Params_t& params) :
			params{ params }
		{}

		template<typename RHS_t>
		State_t operator()(
			const State_t& state,
			const RHS_t& rhs,
			double h) const
		{
			//	using CauchySolver::AbstractSolver;
			State_t k1{ rhs(state, params) };
			State_t k2{ rhs(state + h / 2.0 * k1, params) };
			State_t k3{ rhs(state + h / 2.0 * k2, params) };
			State_t k4{ rhs(state + h * k3, params) };

			return state + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
		}
	};

	template
		<
		typename State_t,
		typename Params_t
		>
		struct RK2
	{
		Params_t params;
		RK2(const Params_t& params) :
			params{ params }
		{}

		template<typename RHS_t>
		State_t operator()(
			const State_t& state,
			const RHS_t& rhs,
			double h) const
		{
			//	using CauchySolver::AbstractSolver;
			State_t k1{ rhs(state, params) };
			State_t k2{ rhs(state + h * k1, params) };

			return state + h / 2.0 * (k1 + k2 );
		}
	};

	struct ConcreteVerifier
	{
		ConcreteVerifier(
			double C, double L,
			double q0, double I0) :
			C{C}, L{L}
		{
			E = (C*q0*q0 + L*I0*I0) / 2.0;
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

