#pragma once
#include <array>

namespace CauchySolver
{
	namespace AbstractSolver
	{
		template<size_t problem_size_t>
		struct State
		{
			template<typename... Args>
			State(Args... state) :
				state{ state }
			{}

			State(const State& state) = default;

			State() = delete;

			auto operator[](size_t i) const
			{
				return state[i];
			}

			std::array<double, problem_size_t>
				state;
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


		template<
			size_t probem_size_t,
			template<size_t problem_size>
		typename State_t,
			typename Params_t>
			struct RK4
		{
			Params_t params;
			RK4(const Params_t& params) :
				params{ params }
			{}

			template<typename Params_t>
			typename RHS_t
				State_t<probem_size_t> operator()(
					const State_t<probem_size_t>& state, 
					const RHS_t<Params_t>& rhs, 
					double h)
			{
				State_t<probem_size_t> k1{ rhs(state, params) };
				State_t<probem_size_t> k2{ rhs(state + h / 2.0 * k1, params) };
				State_t<probem_size_t> k3{ rhs(state + h / 2.0 * k2, params) };
				State_t<probem_size_t> k4{ rhs(state + h * k3, params) };

				return state + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
			}

		};
	}

	

	constexpr size_t problem_size = 4;
	using ConcreteState = AbstractSolver::State<problem_size>;

	struct ConcreteParams
	{
		double omega2;
	};

	struct ConsreteRHS
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

	struct Solver
	{
		ConcreteParams params;

		Solver(const ConcreteParams& params) :
			params{params}
		{}

	};

}
