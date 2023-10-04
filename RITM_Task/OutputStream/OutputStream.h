#pragma once

#include <string>
#include <fstream>
#include <filesystem>

namespace CauchySolver
{
	struct Printer
	{
		using fs = std::filesystem;
		Printer(const std::string& path) :
			path{ path }
		{}

		fs::path path;
	};
} // CauchySolver
