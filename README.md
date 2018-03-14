# CPanel - Unstructured Panel Code

This is an unstructured panel code written in C++. Currently, it models incompressible, 3D, steady or unsteady flows about a single lifting body using constant strength source panels, constant strength double panels, and either prescribed wake panels or vortex particle wakes. The geometry is represented by a triangulated, water-tight surface with an identifiable trailing edge from which the wake will propagate.

It uses [CMake](https://cmake.org/) as the build generator and test running tool. It uses [Google Test](https://github.com/google/googletest) as it's testing framework and [Eigen3]() for linear algebra and matrix math operations. Finally it uses Boost for file I/O.

This project started with a thesis by Chris Satterwhite (2015), [*Development of CPANEL, an Unstructured Panel Code, Using a Modified TLS Velocity Formulation*](https://doi.org/10.15368/theses.2015.135). Specifically, the overall software architecture was developed along with the 3D, steady, prescribed wake panel implementation. More details can be found:

Chris's work was then used as the starting point for Connor Sousa's thesis (2016), [*Unsteady Panel Code Utilizing a Vortex Particle Wake*](http://www.connorsousa.com/thesis/). The prescribed wake implementation was replaced with a vortex particle wake. In addition, the code was modified to be able to model unsteady aerodynamics via surface motion (rotation, translation, etc.) along with software architecture improvements.