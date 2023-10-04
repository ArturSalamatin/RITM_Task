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

				return state + (h / 2.0) * (k1 + k2);
			}
		};
	} // LocalSolver
} // CauchySolver