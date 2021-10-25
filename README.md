# Models-of-Computational-Physics

These projects are three of the most intresting problems assigned during the course of "Metodi Computazioni per la Fisica". These pieces of code are for all those students that for the first time they have to implement a code to solve some challenging physical problems. Note that only some base knowledge of MATLAB is required.

## Fermi Pasta Ulam (FPU) Problem
In 1955 Fermi, Pasta and Ulam had the idea to use and test, one of the first computer ever built the Maniac I, on a simple physical problem: verify the equipartition theorem. For simplicity FPU considered a 1D chain with N points of mass m, one linked to an other with a spring with constant k. They consider also fixed boundary conditions.  In order to verify the equipartition theorem they introduce a non linear term in the equations of motion. They belived that after a certain time, all the modes start to be active. Instead they found that after a long time, called FPU Superperiod, only the initial excited node was active, and so the equipartition theorem seemed not to be verified.
In this code we will plot the famous image that can be found on the original paper LA-1940.

## Poisson's Equation
Poisson's equation is one of the most studied equation in physics and it's one of the first elliptic partial differential equations that a student encounter in his career. It has a broad utility in physics but it's commonly found for the first time in electrodynamics. Indeed with this code we will put a fix distribution of charges on a 2D space and we will calculate the potential and the electric field. The code is actually the implementation of the iterative Gauss-Seidel's algorithm.

## Ising Model
Ising is actually the most cited author in the scientific literature and every physicist has studied his model. In this code we will implement a Metropolis algorithm in order to create a spin-flip dynamics on a 20x20 squared lattice, that will lead the system to the thermalization state. Increasing the temperatue one can indeed verify the phase transition from the magnetic phase to the paramegnetic one. Magnetization and energy per site are as usual the main quantity of interest in this problem.
