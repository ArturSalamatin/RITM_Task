#pragma once
#include <cassert>
#include <iomanip>

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
				assert(i < size());
				return state[i];
			}

			auto& operator[](size_t i)
			{
				assert(i < size());
				return state[i];
			}

			std::array<double, problem_size_t>
				state;

			constexpr static size_t size()
			{
				return problem_size_t;
			}
		};

		// elementwise math operations 
		// with State == array<>
		template<size_t problem_size_t>
		State<problem_size_t> operator+(
			const State<problem_size_t>& lhs,
			const State<problem_size_t>& rhs)
		{
			assert(lhs.size() == rhs.size());
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

		template<typename stream, size_t problem_size_t>
		stream& operator<<(stream& o, const State<problem_size_t>& state)
		{
			for (size_t id = 0; id < state.size(); ++id)
				o
				<< std::scientific
				<< std::setw(17)
				<< std::left
				<< state[id];
		//	o << '\n';

			return o;
		}
	} // AbstractSolver
} // CauchySolver