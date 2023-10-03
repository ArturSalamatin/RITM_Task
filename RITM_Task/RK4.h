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
			State<problem_size_t>& rhs)
		{
			State<problem_size_t> out;
			for (size_t i = 0; i < out; ++i)
				out[i] = lhs[i] + rhs[i];
			return out;
		}

		template<size_t problem_size_t>
		State<problem_size_t> operator*(
			double lhs,
			State<problem_size_t>& rhs)
		{
			State<problem_size_t> out;
			for (size_t i = 0; i < out; ++i)
				out[i] = lhs * rhs[i];
			return out;
		}
		template<size_t problem_size_t>
		State<problem_size_t> operator*(
			State<problem_size_t>& lhs,
			double rhs)
		{
			return rhs * lhs;
		}


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
					double h)
			{
				State_t k1{ rhs(state, params) };
				State_t k2{ rhs(state + h / 2.0 * k1, params) };
				State_t k3{ rhs(state + h / 2.0 * k2, params) };
				State_t k4{ rhs(state + h * k3, params) };

				return state + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
			}

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
	

	constexpr size_t problem_size = 4;
	using ConcreteState = CauchySolver::AbstractSolver::State<problem_size>;

	struct ConcreteParams
	{
		double omega2;
	};

	struct ConcreteRHS
	{
		ConcreteState operator()(const ConcreteState& y, const ConcreteParams& params) const
		{
			return ConcreteState
			{
				y[1],
				-params.omega2 * y[0]
			};
		}
	};

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
			params{params},
			grid{grid},
			init_state{init_state},
			rhs{rhs}
		{
			its_result.reserve(grid.size());
		}

		template<typename LocalSolver_t>
		void solve(const LocalSolver_t& local_solver)
		{
			its_result[0] = init_state;
			for (size_t i = 1; i < grid.size(); ++i)
				its_result[i] = local_solver(result[i - 1], rhs, grid.step(i));
		}

		const std::vector<State_t>& result() const
		{
			if (its_result.size() == 0)
				throw std::exception("One has to call solve() method to compute the result.");
			return its_result;
		}

	protected:
		std::vector<State_t> its_result;
	};

