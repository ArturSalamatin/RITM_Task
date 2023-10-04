#pragma once

namespace CauchySolver
{
	namespace LocalSolver
	{
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
				State_t k1{ rhs(state, params) };
				State_t k2{ rhs(state + (h / 2.0) * k1, params) };
				State_t k3{ rhs(state + (h / 2.0) * k2, params) };
				State_t k4{ rhs(state + h * k3, params) };

				return state + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
			}
		};
	} // LocalSolver
} // CauchySolver
