// RITM_Task.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>

#include "RK4.h"

int main()
{
    double C{ 1.0 };
    double L{ 1.0 };
    double omega2{ C / L };

    ConcreteParams params{ omega2 };

    double U0{ 1.0 };
    double I0{ 0.0 };
    double q0{ U0 / C };
    ConcreteState init_state{U0/C, I0};

    ConcreteRHS rhs{};
    auto init_differential{ rhs(init_state, params) };

    UniformGrid grid{0.0, 100.0, 0.1};

    ConcreteSolver solver{params, grid, init_state, rhs};

    solver.solve(
        RK4<ConcreteState, ConcreteParams>{params}
    );

    ConcreteVerifier verifier{ C, L, q0, I0 };
    verifier.verify(solver);

    std::cout << "Hello World!\n";
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
