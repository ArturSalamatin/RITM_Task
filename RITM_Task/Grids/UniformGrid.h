#pragma once
#include <vector>
#include <cmath>

struct UniformGrid
{
	explicit UniformGrid(double start, double end, size_t node_count) :
		grid(node_count),
		its_step{ (end - start) / (node_count - 1ull) }
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
	size_t size() const
	{
		return grid.size();
	}

	double its_step;
};