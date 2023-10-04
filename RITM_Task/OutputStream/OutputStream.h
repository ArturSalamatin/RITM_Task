#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <filesystem>

#include "ConcreteSolver.h"

namespace CauchySolver
{
	struct Printer
	{
		Printer(
			const std::string& fname = "output.txt",
			const std::string& path = "") :
			path{ path },
			fname{fname}
		{
			if (path == std::string(""))
				abs_path = std::string{ fname };
			else
				abs_path = std::string{ path + std::string{"\\"} + fname};
			
			out_file.open(abs_path);
			if (!out_file.is_open())
				throw std::exception("The output file was not opened correctly.");

			out_file
				<< std::setw(15)
				<< std::left
				<< "Time, sec"

				<< std::setw(17)
				<< std::left
				<< "Calc_Charge, C"

				<< std::setw(17)
				<< std::left
				<< "Calc_Current, A"

				<< std::setw(17)
				<< std::left
				<< "Calc_Energy, J"

				<< std::setw(17)
				<< std::left
				<< "Real_Charge, C"

				<< std::setw(17)
				<< std::left
				<< "Real_Current, A"

				<< std::setw(17)
				<< std::left
				<< "Real_Energy, J"

				<< std::setw(21)
				<< std::left
				<< "Abs_Energy_diff, J"

				<< std::setw(21)
				<< std::left
				<< "Rel_Energy_diff, -"

				<< '\n';
		}

		void print(const ConcreteVerifier& verifier, const ConcreteSolver& solver)
		{
			for (size_t t_id = 0; t_id < solver.grid.size(); ++t_id)
			{
				out_file

					<< std::scientific
					<< std::setw(15)
					<< std::left
					<< solver.grid.grid[t_id]

					<< solver.result(t_id)

					<< std::scientific
					<< std::setw(17)
					<< std::left
					<< verifier.calc_energy[t_id]

					<< std::scientific
					<< std::setw(17)
					<< std::left
					<< verifier.real_charge[t_id]

					<< std::scientific
					<< std::setw(17)
					<< std::left
					<< verifier.real_current[t_id]

					<< std::scientific
					<< std::setw(17)
					<< std::left
					<< verifier.E

					<< std::scientific
					<< std::setw(21)
					<< std::left
					<< verifier.abs_diff[t_id]

					<< std::scientific
					<< std::setw(21)
					<< std::left
					<< verifier.rel_diff[t_id]

					<< '\n';
			}
		}

		std::string path;
		std::string fname;
		std::string abs_path;

		std::ofstream out_file;
	};
} // CauchySolver
