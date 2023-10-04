#pragma once
#include <array>
#include <exception>

namespace CauchySolver
{
	namespace AbstractSolver
	{
		template<
			typename State_t,
			typename Params_t,
			typename Grid_t,
			typename RHS_t>
			struct AbstractSolver
		{
			Params_t params;
			Grid_t grid;
			State_t init_state;
			RHS_t rhs;

			AbstractSolver(
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